from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import yaml

from wam2layers.preprocessing.preprocessing import (accumulation_to_flux, calculate_humidity, insert_level, interpolate,
                           sortby_ndarray)
from wam2layers.analysis.checks import check_input



def load_data(variable, date, levels=None):
    """Load data for given variable and date."""
    # TODO: remove hardcoded filename pattern
    filepath = Path(config["input_folder"]) / f"FloodCase_202107_{variable}.nc"
    da = xr.open_dataset(filepath)[variable]

    # Include midnight of the next day (if available)
    extra = date + pd.Timedelta(days=1)

    if levels is not None:
        return da.sel(time=slice(date, extra)).sel(level=levels)
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


def get_dp_pressurelevels(q, u, v, ps, qs, us, vs):
    """Get dp with consistent u, v, q for ERA5 pressure level data."""
    p = u.level.broadcast_like(u) * 100  # Pa

    # Insert top of atmosphere values
    # Assume wind at top same as values at lowest pressure, humidity at top 0
    u = insert_level(u, u.isel(level=0), 0)
    v = insert_level(v, v.isel(level=0), 0)
    q = insert_level(q, 0, 0)
    p = insert_level(p, 0, 0)

    # Insert surface level values (at a high dummy pressure value)
    u = insert_level(u, us, 110000)
    v = insert_level(v, vs, 110000)
    q = insert_level(q, qs, 110000)
    p = insert_level(p, ps, 110000)

    # Sort arrays by pressure (ascending)
    u.values = sortby_ndarray(u.values, p.values, axis=1)
    v.values = sortby_ndarray(v.values, p.values, axis=1)
    q.values = sortby_ndarray(q.values, p.values, axis=1)
    p.values = sortby_ndarray(p.values, p.values, axis=1)

    # Insert boundary level values (at a ridiculous dummy pressure value)
    p_boundary = 0.72878581 * np.array(ps) + 7438.803223
    u = insert_level(u, interpolate(p_boundary, p, u), 150000)
    v = insert_level(v, interpolate(p_boundary, p, v), 150000)
    q = insert_level(q, interpolate(p_boundary, p, q), 150000)
    p = insert_level(p, p_boundary, 150000)

    # Sort arrays by pressure once more (ascending)
    u.values = sortby_ndarray(u.values, p.values, axis=1)
    v.values = sortby_ndarray(v.values, p.values, axis=1)
    q.values = sortby_ndarray(q.values, p.values, axis=1)
    p.values = sortby_ndarray(p.values, p.values, axis=1)

    # Reset level coordinate as its values have become meaningless
    nlev = u.level.size
    u = u.assign_coords(level=np.arange(nlev))
    v = v.assign_coords(level=np.arange(nlev))
    q = q.assign_coords(level=np.arange(nlev))
    p = p.assign_coords(level=np.arange(nlev))

    # Calculate pressure jump
    dp = p.diff("level")
    assert np.all(dp > 0), "Pressure levels should increase monotonically"

    # Interpolate to midpoints
    midpoints = 0.5 * (u.level.values[1:] + u.level.values[:-1])
    dp = dp.assign_coords(level=midpoints)
    u = u.interp(level=midpoints)
    v = v.interp(level=midpoints)
    q = q.interp(level=midpoints)
    p = p.interp(level=midpoints)

    # mask values below surface
    above_surface = p < np.array(ps)[:, None, :, :]
    u = u.where(above_surface)
    v = v.where(above_surface)
    q = q.where(above_surface)
    p = p.where(above_surface)

    return dp, p, q, u, v, p_boundary


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


level_type = "pressure_levels"
config = parse_config("../../cases/era5_2021_local.yaml")
for date in config["datelist"]:
    print(date)
    precip, evap = preprocess_precip_and_evap(date)

    # 4d fields
    q = load_data("q", date)  # in kg kg-1
    u = load_data("u", date)  # in m/s
    v = load_data("v", date)  # in m/s
    sp = load_data("sp", date)  # in Pa

    if level_type == "pressure_levels":
        d2m = load_data("d2m", date)  # Dew point in K
        q2m = calculate_humidity(d2m, sp)  # kg kg-1
        u10 = load_data("u10", date)  # in m/s
        v10 = load_data("v10", date)  # in m/s
        dp, p, q, u, v, pb = get_dp_pressurelevels(q, u, v, sp, q2m, u10, v10)

    g = 9.80665  # gravitational accelleration [m/s2]
    cwv = q * dp / g  # column water vapor (kg/m2)

    if config["vertical_integral_available"]:
        # calculate column water instead of column water vapour
        tcw = load_data("tcw", date)  # kg/m2
        cw = (tcw / cwv.sum(dim="level")) * cwv  # column water (kg/m2)
        # TODO: warning if cw >> cwv
    else:
        # fluxes will be calculated based on the column water vapour
        cw = cwv

    # Integrate fluxes and states to upper and lower layer
    if level_type == "pressure_levels":
        upper_layer = p < pb[:, None, :, :]
        lower_layer = pb[:, None, :, :] < p

    # Vertically integrate state over two layers
    s_lower = cw.where(lower_layer).sum(dim="level")
    s_upper = cw.where(upper_layer).sum(dim="level")

    # Determine the fluxes
    fx = u * cw  # eastward atmospheric moisture flux (kg m-1 s-1)
    fy = v * cw  # northward atmospheric moisture flux (kg m-1 s-1)

    # Vertically integrate fluxes over two layers
    fx_lower = fx.where(lower_layer).sum(dim="level")  # kg m-1 s-1
    fy_lower = fy.where(lower_layer).sum(dim="level")  # kg m-1 s-1
    fx_upper = fx.where(upper_layer).sum(dim="level")  # kg m-1 s-1
    fy_upper = fy.where(upper_layer).sum(dim="level")  # kg m-1 s-1

    # Combine everything into one dataset
    d2m = xr.Dataset(
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
    check_input(d2m)

    # Save preprocessed data
    filename = f"{date.strftime('%Y-%m-%d')}_fluxes_storages.nc"
    output_path = config["output_dir"] / filename
    d2m.to_netcdf(output_path)
