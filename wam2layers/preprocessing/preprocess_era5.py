from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import yaml

from preprocessing import (get_grid_info, insert_level, interpolate,
                           sortby_ndarray, calculate_humidity)

# Set constants
g = 9.80665  # [m/s2]
density_water = 1000  # [kg/m3]

# Read case configuration
with open("../../cases/era5_2013.yaml") as f:
    config = yaml.safe_load(f)

# Create the preprocessed data folder if it does not exist yet
output_dir = Path(config["preprocessed_data_folder"]).expanduser()
output_dir.mkdir(exist_ok=True, parents=True)


def load_data(variable, date):
    """Load data for given variable and date."""
    filepath = Path(config["input_folder"]) / f"FloodCase_201305_{variable}.nc"
    da = xr.open_dataset(filepath)[variable]

    # Include midnight of the next day (if available)
    extra = date + pd.Timedelta(days=1)
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
    q = load_data("q", date)  # in kg kg-1
    u = load_data("u", date)  # in m/s
    v = load_data("v", date)  # in m/s

    # Precipitation and evaporation
    evap = load_data("e", date)  # in m (accumulated hourly)
    cp = load_data("cp", date)  # convective precipitation in m (accumulated hourly)
    lsp = load_data("lsp", date)  # large scale precipitation in m (accumulated hourly)
    precip = cp + lsp

    # TODO: not used
    tcw = load_data("tcw", date)  # kg/m2
    
    # surface variables
    p_surf = load_data("sp", date)  # in Pa
    d_surf = load_data("d2m", date)  # Dew point in K
    u_surf = load_data("u10", date)  # in m/s
    v_surf = load_data("v10", date)  # in m/s
    q_surf = calculate_humidity(d_surf, p_surf)  # kg kg-1

    # Get grid info
    time = u.time.values
    lat = u.latitude.values
    lon = u.longitude.values
    a_gridcell, l_ew_gridcell, l_mid_gridcell = get_grid_info(u)

    # Calculate volumes
    # TODO: shall we do this in the tracking?
    evap *= a_gridcell[None, :, None]  # m3
    precip *= a_gridcell[None, :, None]  # m3

    # Transfer negative (originally positive) values of evap to precip
    precip = np.maximum(precip, 0) + np.maximum(evap, 0)
    # Change sign convention to all positive,
    evap = np.abs(np.minimum(evap, 0))

    # Create pressure array with the same dimensions as u, q, and v and convert to Pascal from hPa
    p = u.level.broadcast_like(u) * 100 # Pa

    # Insert top of atmosphere values
    u = insert_level(u, u.isel(level=0), 0) # assume wind at top same as value at lowest pressure
    v = insert_level(v, v.isel(level=0), 0) # assume wind at top same as value at lowest pressure
    q = insert_level(q, 0, 0) # assume pressure at top 0
    p = insert_level(p, 0, 0)

    # Insert surface level values (at a high dummy pressure value)
    u = insert_level(u, u_surf, 110000)
    v = insert_level(v, v_surf, 110000)
    q = insert_level(q, q_surf, 110000)
    p = insert_level(p, p_surf, 110000)
    
    # set pressure value below surface to surface pressure
    below_surface = p > np.array(p_surf)[:, None, :, :]
    # TODO: insert a line that replaces all pressure values where the pressure is below surface with surface pressure (without a loop), something like this, but does not work (array dimension mismatch):
    #p[p.where(below_surface)] = p_surf

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
    dp.assign_coords(level=midpoints)
    
    # water vapor voxels #TODO: ERROR the lines below results for unknown reason in cwv with shape=(25, 0, 121, 321)
    cwv = q * dp / g  #  # column water vapor (kg/m2)
    
    if config["vertical_integral_available"] == True:
        # calculate column water instead of column water vapour
        tcw = load_data("tcw", date)  # kg/m2
        cw = (tcw / cwv.sum(dim="lev")) * cwv  # column water (kg/m2)
    else:
        # calculate the fluxes based on the column water vapour
        cw = cwv
    
    # Integrate fluxes and states to upper and lower layer
    upper_layer = p < p_boundary[:, None, :, :]
    lower_layer = ~upper_layer
    lower_layer = (dp.level < sp / 100) & (dp.level > P_boundary / 100)
             
    # Vertically integrate state over two layers
    s_lower = (
        cw.where(lower_layer).sum(dim="lev")
        * a_gridcell[None, :, None]
        / density_water
    )  # m3
    s_upper = (
        cw.where(upper_layer).sum(dim="lev")
        * a_gridcell[None, :, None]
        / density_water
    )  # m3

    # Determine the fluxes
    fx = u * cw  # eastward atmospheric moisture flux (kg m-1 s-1)
    fy = v * cw  # northward atmospheric moisture flux (kg m-1 s-1)

    # Vertically integrate fluxes over two layers
    fx_lower = fx.where(lower_layer).sum(dim="lev")  # kg m-1 s-1
    fy_lower = fy.where(lower_layer).sum(dim="lev")  # kg m-1 s-1

    fx_upper = fx.where(upper_layer).sum(dim="lev")  # kg m-1 s-1
    fy_upper = fy.where(upper_layer).sum(dim="lev")  # kg m-1 s-1

    # Save preprocessed data
    filename = f"{date.strftime('%Y-%m-%d')}_fluxes_storages.nc"
    output_path = output_dir / filename

    xr.Dataset(
        {  # TODO: would be nice to add coordinates and units as well
            "fx_upper": fx_upper,
            "fy_upper": fy_upper,
            "fx_lower": fx_lower,
            "fy_lower": fy_lower,
            "s_upper": s_upper,
            "s_lower": s_lower,
            "evap": evap,
            "precip": precip,
        }
    ).to_netcdf(output_path)
