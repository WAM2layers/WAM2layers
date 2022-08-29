from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import yaml

from .preprocessing import (calculate_humidity, insert_level, interpolate,
                           sortby_ndarray)
from wam2layers.analysis.checks import check_input

# Set constants
g = 9.80665  # [m/s2]

# Read case configuration
with open("../../cases/era5_2021.yaml") as f:
    config = yaml.safe_load(f)

# Create the preprocessed data folder if it does not exist yet
output_dir = Path(config["preprocessed_data_folder"]).expanduser()
output_dir.mkdir(exist_ok=True, parents=True)


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


datelist = pd.date_range(
    start=config["preprocess_start_date"],
    end=config["preprocess_end_date"],
    freq="d",
    inclusive="left",
)

for date in datelist[:]:
    print(date)

    # 4d fields
    u = load_data("u", date)  # in m/s
    v = load_data("v", date)  # in m/s
    q = load_data("q", date)  # in kg kg-1

    # Precipitation and evaporation
    evap = load_data("e", date)  # in m (accumulated hourly)
    cp = load_data("cp", date)  # convective precip in m (accumulated hourly)
    lsp = load_data("lsp", date)  # large scale precip in m (accumulated)
    precip = cp + lsp

    # 3d fields
    p_surf = load_data("sp", date)  # in Pa
    d_surf = load_data("d2m", date)  # Dew point in K
    u_surf = load_data("u10", date)  # in m/s
    v_surf = load_data("v10", date)  # in m/s
    q_surf = calculate_humidity(d_surf, p_surf)  # kg kg-1

    # Transfer negative (originally positive) values of evap to precip
    precip = np.maximum(precip, 0) + np.maximum(evap, 0)
    # Change sign convention to all positive,
    evap = np.abs(np.minimum(evap, 0))

    # Convert precip and evap from accumulations to fluxes
    density = 1000  # [kg/m3]
    timestep = pd.Timedelta(precip.time.diff('time')[0].values)
    nseconds = timestep.total_seconds()
    precip = (density * precip / nseconds).assign_attrs(units="kg m-2 s-1")
    evap = (density * evap / nseconds).assign_attrs(units="kg m-2 s-1")

    # Valid time are now midway throught the timestep, better to be consistent
    # Extrapolation introduces a small inconsistency at the last midnight...
    original_time = precip.time
    midpoint_time = original_time - timestep / 2
    precip["time"] = precip.time - timestep / 2
    evap["time"] = precip.time - timestep / 2
    precip.interp(time=original_time, kwargs={"fill_value": "extrapolate"})
    evap.interp(time=original_time, kwargs={"fill_value": "extrapolate"})


    # Create pressure array with the same dimensions as u, q, and v and convert
    # from hPa to Pa
    p = u.level.broadcast_like(u) * 100  # Pa

    # Insert top of atmosphere values
    # Assume wind at top same as values at lowest pressure, humidity at top 0
    u = insert_level(u, u.isel(level=0), 0)
    v = insert_level(v, v.isel(level=0), 0)
    q = insert_level(q, 0, 0)
    p = insert_level(p, 0, 0)

    # Insert surface level values (at a high dummy pressure value)
    u = insert_level(u, u_surf, 110000)
    v = insert_level(v, v_surf, 110000)
    q = insert_level(q, q_surf, 110000)
    p = insert_level(p, p_surf, 110000)

    # Sort arrays by pressure (ascending)
    u.values = sortby_ndarray(u.values, p.values, axis=1)
    v.values = sortby_ndarray(v.values, p.values, axis=1)
    q.values = sortby_ndarray(q.values, p.values, axis=1)
    p.values = sortby_ndarray(p.values, p.values, axis=1)

    # Insert boundary level values (at a ridiculous dummy pressure value)
    p_boundary = 0.72878581 * np.array(p_surf) + 7438.803223
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
    u = u.interp(level=midpoints)
    v = v.interp(level=midpoints)
    q = q.interp(level=midpoints)
    p = p.interp(level=midpoints)
    dp = dp.assign_coords(level=midpoints)

    # mask values below surface
    above_surface = p < np.array(p_surf)[:, None, :, :]
    u = u.where(above_surface)
    v = v.where(above_surface)
    q = q.where(above_surface)
    p = p.where(above_surface)

    # water vapor voxels
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
    upper_layer = p < p_boundary[:, None, :, :]
    lower_layer = p_boundary[:, None, :, :] < p

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
    output_path = output_dir / filename
    ds.to_netcdf(output_path)
