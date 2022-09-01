from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import yaml

from wam2layers.analysis.checks import check_input
from wam2layers.preprocessing.preprocessing import accumulation_to_flux


def load_data(variable, date, levels=None):
    """Load data for given variable and date."""
    # TODO: remove hardcoded filename pattern
    # TODO: make sure lev coordinate is named consistently
    prefix = "FloodCase_202107_ml" if levels is not None else "FloodCase_202107"
    filepath = Path(config["input_folder"]) / f"{prefix}_{variable}.nc"
    da = xr.open_dataset(filepath)[variable]

    # Include midnight of the next day (if available)
    extra = date + pd.Timedelta(days=1)

    if levels is not None:
        return da.sel(time=slice(date, extra)).sel(lev=levels)
    return da.sel(time=slice(date, extra))


def preprocess_precip_and_evap(date):
    """Load and pre-process precipitation and evaporation."""
    # All incoming units are accumulations (in m) since previous time step
    evap = load_data("e", date)
    cp = load_data("cp", date)  # convective precipitation
    lsp = load_data("lsp", date)  # large scale precipitation
    precip = (cp + lsp)

    # Transfer negative (originally positive) values of evap to precip
    precip = np.maximum(precip, 0) + np.maximum(evap, 0)
    precip = precip

    # Change sign convention to all positive,
    evap = np.abs(np.minimum(evap, 0))

    precip = accumulation_to_flux(precip)
    evap = accumulation_to_flux(evap)

    return precip, evap


def get_edges(era5_modellevels):
    """Get the values of a and b at the edges of a subset of ERA5 modellevels."""
    # Load a and b coefficients
    df = pd.read_csv("tableERA5model_to_pressure.csv")
    a = df["a [Pa]"].to_xarray().rename(index="lev")
    b = df["b"].to_xarray().rename(index="lev")

    # Calculate a and b at mid levels (model levels)
    a_full = ((a[1:] + a[:-1].values) / 2.0).sel(lev=era5_modellevels)
    b_full = ((b[1:] + b[:-1].values) / 2.0).sel(lev=era5_modellevels)

    # Interpolate to get parameters at edges between selected levels
    a_edge = xr.concat([a[0], (a_full[1:].values + a_full[:-1]) / 2, a[-1]], dim="lev")
    b_edge = xr.concat([b[0], (b_full[1:].values + b_full[:-1]) / 2, b[-1]], dim="lev")

    return a_edge, b_edge


def get_dp_modellevels(sp, modellevels):
    """Calculate pressure jump over subset of era5 modellevels."""
    sp_broadcast = sp.expand_dims({"lev": modellevels}, axis=1)
    a_edge, b_edge = get_edges(modellevels)
    p_edges = a_edge + b_edge * sp_broadcast  # in Pa

    # calculate the difference between the pressure levels
    dp_modellevels = p_edges.diff(dim="lev")  # in Pa

    return dp_modellevels


def parse_config(config_file):
    """Read and parse case configuration file."""
    with open(config_file) as f:
        config = yaml.safe_load(f)

    # Create the preprocessed data folder if it does not exist yet
    config["output_dir"] = Path(config["preprocessed_data_folder"]).expanduser()
    config["output_dir"].mkdir(exist_ok=True, parents=True)

    config["datelist"] = pd.date_range(
        start=config["preprocess_start_date"],
        end=config["preprocess_end_date"],
        freq="d",
        inclusive="left",
    )
    return config


level_type = "model_levels"
config = parse_config("../../cases/era5_2021_local.yaml")
for date in config["datelist"]:
    print(date)
    precip, evap = preprocess_precip_and_evap(date)

    # 4d fields
    modellevels = config["modellevels"]
    q = load_data("q", date, levels=modellevels)  # in kg kg-1
    u = load_data("u", date, levels=modellevels)  # in m/s
    v = load_data("v", date, levels=modellevels)  # in m/s
    sp = load_data("sp", date)  # in Pa

    if level_type == "model_levels":
        dp = get_dp_modellevels(sp, modellevels)


    # Calculate column water vapour
    g = 9.80665  # gravitational accelleration [m/s2]
    cwv = q * dp / g  # (kg/m2)

    if config["vertical_integral_available"]:
        # calculate column water instead of column water vapour
        tcw = load_data("tcw", date)  # kg/m2
        cw = (tcw / cwv.sum(dim="lev")) * cwv  # column water (kg/m2)
        # TODO: warning if cw >> cwv
    else:
        # fluxes will be calculated based on the column water vapour
        cw = cwv

    # Integrate fluxes and states to upper and lower layer
    # TODO: Check if this is a reasonable choice for boundary
    if level_type == "model_levels":
        boundary = 111
        lower_layer = dp.lev > boundary
        upper_layer = ~lower_layer

    # Vertically integrate state over two layers
    s_lower = cw.where(lower_layer).sum(dim="lev")
    s_upper = cw.where(upper_layer).sum(dim="lev")

    # Determine the fluxes
    fx = (u * cw)  # eastward atmospheric moisture flux (kg m-1 s-1)
    fy = (v * cw)  # northward atmospheric moisture flux (kg m-1 s-1)

    # Vertically integrate fluxes over two layers
    fx_lower = fx.where(lower_layer).sum(dim="lev")  # kg m-1 s-1
    fy_lower = fy.where(lower_layer).sum(dim="lev")  # kg m-1 s-1
    fx_upper = fx.where(upper_layer).sum(dim="lev")  # kg m-1 s-1
    fy_upper = fy.where(upper_layer).sum(dim="lev")  # kg m-1 s-1

    # Combine everything into one dataset
    ds = xr.Dataset(
        {  # TODO: would be nice to add coordinates and units as well
            "fx_upper": fx_upper.assign_attrs(units="kg m-1 s-1"),
            "fy_upper": fy_upper.assign_attrs(units="kg m-1 s-1"),
            "fx_lower": fx_lower.assign_attrs(units="kg m-1 s-1"),
            "fy_lower": fy_lower.assign_attrs(units="kg m-1 s-1"),
            "s_upper": s_upper.assign_attrs(units="kg m-2"),
            "s_lower": s_lower.assign_attrs(units="kg m-2"),
            "evap": evap,
            "precip": precip,
        }
    )

    # Verify that the data meets all the requirements for the model
    check_input(ds)

    # Save preprocessed data
    filename = f"{date.strftime('%Y-%m-%d')}_fluxes_storages.nc"
    output_path = config["output_dir"] / filename
    ds.to_netcdf(output_path)
