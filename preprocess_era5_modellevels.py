import os

import pandas as pd
import xarray as xr
import yaml
import numpy as np
from pathlib import Path

from preprocessing import get_grid_info, get_stable_fluxes, get_vertical_transport


# Read case configuration
with open("cases/era5_2021.yaml") as f:
    config = yaml.safe_load(f)

# Create the preprocessed data folder if it does not exist yet
output_dir = Path(config["preprocessed_data_folder"]).expanduser()
output_dir.mkdir(exist_ok=True, parents=True)

# Set constants
g = 9.80665  # [m/s2]
density_water = 1000  # [kg/m3]

# Preselection of model levels
modellevels = config["modellevels"]
print("Number of model levels:", len(modellevels))

# Load a and b coefficients
df = pd.read_csv("tableERA5model_to_pressure.csv")
a = df["a [Pa]"].to_xarray().rename(index="lev")  # .sel(lev=modellevels)
b = df["b"].to_xarray().rename(index="lev")  # .sel(lev=modellevels)

# Calculate a and b at mid levels (model levels)
a_full = ((a[1:] + a[:-1].values) / 2.0).sel(lev=modellevels)
b_full = ((b[1:] + b[:-1].values) / 2.0).sel(lev=modellevels)

# Construct a and b at edges for selected levels
a_edge = xr.concat([a[0], (a_full[1:].values + a_full[:-1]) / 2.0, a[-1]], dim="lev")
b_edge = xr.concat([b[0], (b_full[1:].values + b_full[:-1]) / 2.0, b[-1]], dim="lev")


def load_surface_data(variable, date):
    """Load data for given variable and date."""
    filename = f"FloodCase_202107_{variable}.nc"
    filepath = os.path.join(config["input_folder"], filename)
    da = xr.open_dataset(filepath)[variable]

    # Include midnight of the next day (if available)
    extra = date + pd.Timedelta(days=1)
    return da.sel(time=slice(date, extra))


def load_modellevel_data(variable, date):
    """Load model level data for given variable and date."""
    filename = f"FloodCase_202107_ml_{variable}.nc"
    filepath = os.path.join(config["input_folder"], filename)
    da = xr.open_dataset(filepath)[variable]

    # Include midnight of the next day (if available)
    extra = date + pd.Timedelta(days=1)
    return da.sel(time=slice(date, extra)).sel(lev=modellevels)


datelist = pd.date_range(
    start=config["preprocess_start_date"],
    end=config["preprocess_end_date"],
    freq="d",
    inclusive="left",
)

for date in datelist[:]:
    print(date)

    # Load data
    u = load_modellevel_data("u", date)
    v = load_modellevel_data("v", date)
    q = load_modellevel_data("q", date)
    sp = load_surface_data("sp", date)  # in Pa
    evap = load_surface_data("e", date)
    cp = load_surface_data("cp", date)
    lsp = load_surface_data("lsp", date)
    precip = cp + lsp

    # Get grid info
    lat = u.latitude.values
    lon = u.longitude.values
    a_gridcell, l_ew_gridcell, l_mid_gridcell = get_grid_info(lat, lon)

    # Calculate volumes
    evap *= a_gridcell  # m3
    precip *= a_gridcell  # m3

    # Transfer negative (originally positive) values of evap to precip
    precip = np.maximum(precip, 0) + np.maximum(evap, 0)
    # Change sign convention to all positive,
    evap = np.abs(np.minimum(evap, 0))

    # Convert model levels to pressure values
    sp_modellevels2 = sp.expand_dims({"lev": modellevels}, axis=1)
    p_modellevels = a_edge + b_edge * sp_modellevels2  # in Pa

    # calculate the difference between the pressure levels
    dp_modellevels = p_modellevels.diff(dim="lev")  # in Pa

    # Split in 2 layers
    # To do: Check if this is a reasonable choice
    boundary = 111  # set boundary model level 111

    # Integrate fluxes and states to upper and lower layer
    lower_layer = dp_modellevels.lev > boundary
    upper_layer = ~lower_layer

    if config["vertical_integral_available"] == True:
        # load total water column
        tcw = load_surface_data("tcw", date)  # kg/m2

        # Determine the states
        cwv = q * dp_modellevels / g  #  # column water vapor (kg/m2)

        # calculate column water instead of column water vapour
        cw = (tcw / cwv.sum(dim="lev")) * cwv  # column water (kg/m2)

        # Determine the fluxes
        fx = u * cw  # eastward atmospheric moisture flux (kg m-1 s-1)
        fy = v * cw  # northward atmospheric moisture flux (kg m-1 s-1)

        # no correction for fluxes as it is hard to decide how to distribute the correction over the two layers
        # no correction for fluxes is applied but we use the column water to calculate them

        # Vertically integrate state over two layers
        s_lower = (
            cw.where(lower_layer).sum(dim="lev")
            * a_gridcell[np.newaxis, :]
            / density_water
        )  # m3
        s_upper = (
            cw.where(upper_layer).sum(dim="lev")
            * a_gridcell[np.newaxis, :]
            / density_water
        )  # m3

    else:  # calculate the fluxes based on the column water vapour

        # Determine the states
        cwv = q * dp_modellevels / g  #  # column water vapor (kg/m2)

        # Determine the fluxes
        fx = u * cwv  # eastward atmospheric moisture flux (kg m-1 s-1)
        fy = v * cwv  # northward atmospheric moisture flux (kg m-1 s-1)

        # Vertically integrate state over two layers

        s_lower = (
            cwv.where(lower_layer).sum(dim="lev")
            * a_gridcell[np.newaxis, :]
            / density_water
        )  # m3
        s_upper = (
            cwv.where(upper_layer).sum(dim="lev")
            * a_gridcell[np.newaxis, :]
            / density_water
        )  # m3

    # Vertically integrate fluxes over two layers
    fx_lower = fx.where(lower_layer).sum(dim="lev")  # kg m-1 s-1
    fy_lower = fy.where(lower_layer).sum(dim="lev")  # kg m-1 s-1

    fx_upper = fx.where(upper_layer).sum(dim="lev")  # kg m-1 s-1
    fy_upper = fy.where(upper_layer).sum(dim="lev")  # kg m-1 s-1

    # Change units to m3, based on target frequency (not incoming frequency!)
    target_freq = config["target_frequency"]
    total_seconds = pd.Timedelta(target_freq).total_seconds()
    fx_upper *= total_seconds * (l_ew_gridcell / density_water)
    fx_lower *= total_seconds * (l_ew_gridcell / density_water)
    fy_upper *= total_seconds * (l_mid_gridcell[None, :, None] / density_water)
    fy_lower *= total_seconds * (l_mid_gridcell[None, :, None] / density_water)

    # Put data on a smaller time step...
    time = s_upper.time.values
    newtime = pd.date_range(time[0], time[-1], freq="15Min")[:-1]
    s_upper = s_upper.interp(time=newtime).values
    s_lower = s_lower.interp(time=newtime).values

    # ... fluxes on the edges instead of midpoints
    newtime = newtime[:-1] + pd.Timedelta("6Min") / 2
    fx_upper = fx_upper.interp(time=newtime).values
    fy_upper = fy_upper.interp(time=newtime).values

    fx_lower = fx_lower.interp(time=newtime).values
    fy_lower = fy_lower.interp(time=newtime).values

    precip = (precip.reindex(time=newtime, method="bfill") / 4).values
    evap = (evap.reindex(time=newtime, method="bfill") / 4).values

    # Stabilize horizontal fluxes
    fx_upper, fy_upper = get_stable_fluxes(fx_upper, fy_upper, s_upper)
    fx_lower, fy_lower = get_stable_fluxes(fx_lower, fy_lower, s_lower)

    # Determine the vertical moisture flux
    f_vert = get_vertical_transport(
        fx_upper,
        fx_lower,
        fy_upper,
        fy_lower,
        evap,
        precip,
        s_upper,
        s_lower,
        config["periodic_boundary"],
        config["kvf"],
    )

    # Save preprocessed data
    # Note: fluxes (dim: time) are at the edges of the timesteps,
    # while states (dim: time2) are at the midpoints and include next midnight
    # so the first state from day 2 will overlap with the last flux from day 1
    filename = f"{date.strftime('%Y-%m-%d')}_fluxes_storages.nc"
    output_path = os.path.join(config["preprocessed_data_folder"], filename)
    xr.Dataset(
        {  # TODO: would be nice to add coordinates and units as well
            "fx_upper": (["time", "lat", "lon"], fx_upper),
            "fy_upper": (["time", "lat", "lon"], fy_upper),
            "fx_lower": (["time", "lat", "lon"], fx_lower),
            "fy_lower": (["time", "lat", "lon"], fy_lower),
            "s_upper": (["time2", "lat", "lon"], s_upper),
            "s_lower": (["time2", "lat", "lon"], s_lower),
            "f_vert": (["time", "lat", "lon"], f_vert),
            "evap": (["time", "lat", "lon"], evap),
            "precip": (["time", "lat", "lon"], precip),
        }
    ).to_netcdf(output_path)
