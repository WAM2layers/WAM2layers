import os
import numpy as np
import pandas as pd
import numpy as np
import xarray as xr
import yaml


from preprocessing import (get_grid_info, get_stable_fluxes,
                           get_vertical_transport)


# Set constants
g = 9.80665  # [m/s2]
density_water = 1000  # [kg/m3]

# Read case configuration
with open("cases/era5_2013.yaml") as f:
    config = yaml.safe_load(f)


def load_data(variable, date):
    """Load data for given variable and date."""
    filename = f"FloodCase_201305_{variable}.nc"
    filepath = os.path.join(config["input_folder"], filename)
    da = xr.open_dataset(filepath)[variable]

    # Include midnight of the next day (if available)
    extra = date + pd.Timedelta(days=1)
    return da.sel(time=slice(date, extra))


datelist = pd.date_range(
    start=config["preprocess_start_date"], end=config["preprocess_end_date"], freq="d", inclusive="left"
)

for date in datelist[:]:
    print(date)

    # Load data
    u = load_data("u", date) #in m/s
    v = load_data("v", date) #in m/s
    q = load_data("q", date) #in kg kg-1
    sp = load_data("sp", date) #in Pa
    evap = load_data("e", date) #in m (accumulated hourly)
    cp = load_data("cp", date) #convective precipitation in m (accumulated hourly)
    lsp = load_data("lsp", date) #large scale precipitation in m (accumulated hourly)
    precip = cp + lsp
    tcw = load_data("tcw", date) # kg/m2

    # Get grid info
    lat = u.latitude.values
    lon = u.longitude.values
    a_gridcell, l_es_gridcell, l_mid_gridcell = get_grid_info(lat, lon)

    # Calculate volumes
    evap *= a_gridcell  # m3
    precip *= a_gridcell  # m3

    # Transfer negative (originally positive) values of evap to precip
    precip = np.maximum(precip, 0) + np.maximum(evap, 0)
    # Change sign convention to all positive,
    evap = np.abs(np.minimum(evap, 0))

    # Create pressure array
    levels = q.level #in hPa
    p = levels.broadcast_like(u) * 100  # Pa

    # Interpolate to new levels
    edges = 0.5 * (levels.values[1:] + levels.values[:-1])
    u = u.interp(level=edges)
    v = v.interp(level=edges)
    q = q.interp(level=edges)

    # Calculate pressure jump
    dp = p.diff(dim="level")
    dp["level"] = edges

    # Determine the fluxes and states
    fx = u * q * dp / g  # eastward atmospheric moisture flux (kg*m-1*s-1)
    fy = v * q * dp / g  # northward atmospheric moisture flux (#kg*m-1*s-1)
    cwv = q * dp / g * a_gridcell / density_water  # column water vapor (m3)

    # Split in 2 layers
    P_boundary = 0.72878581 * sp + 7438.803223
    lower_layer = (dp.level < sp / 100) & (dp.level > P_boundary / 100)
    upper_layer = dp.level < P_boundary / 100

    # Integrate fluxes and state
    fx_lower = fx.where(lower_layer).sum(dim="level") #kg*m-1*s-1
    fy_lower = fy.where(lower_layer).sum(dim="level") #kg*m-1*s-1
    s_lower = cwv.where(lower_layer).sum(dim="level") #m3

    fx_upper = fx.where(upper_layer).sum(dim="level") #kg*m-1*s-1
    fy_upper = fy.where(upper_layer).sum(dim="level") #kg*m-1*s-1
    s_upper = cwv.where(upper_layer).sum(dim="level") #m3

    print(
        "Check calculation water vapor, this value should be zero:",
        (cwv.sum(dim="level") - (s_upper + s_lower)).sum().values,
    )

    tcwm3 = tcw * a_gridcell[np.newaxis, :] / density_water  # m3

    # Change units to m3, based on target frequency (not incoming frequency!)
    target_freq = config['target_frequency']
    total_seconds = pd.Timedelta(target_freq).total_seconds()
    fx_upper *= total_seconds * (l_es_gridcell / density_water)
    fx_lower *= total_seconds * (l_es_gridcell / density_water)
    fy_upper *= total_seconds * (l_mid_gridcell[None, :, None] / density_water)
    fy_lower *= total_seconds * (l_mid_gridcell[None, :, None] / density_water)

    # Increased time resolution; states at midpoints, fluxes at the edges
    old_time = s_upper.time.values
    newtime_states = pd.date_range(old_time[0], old_time[-1], freq=target_freq)
    newtime_fluxes = newtime_states[:-1] + pd.Timedelta(target_freq) / 2

    s_upper = s_upper.interp(time=newtime_states).values
    s_lower = s_lower.interp(time=newtime_states).values
    fx_upper = fx_upper.interp(time=newtime_fluxes).values
    fy_upper = fy_upper.interp(time=newtime_fluxes).values
    fx_lower = fx_lower.interp(time=newtime_fluxes).values
    fy_lower = fy_lower.interp(time=newtime_fluxes).values
    precip = (precip.reindex(time=newtime_fluxes, method="bfill") / 4).values
    evap = (evap.reindex(time=newtime_fluxes, method="bfill") / 4).values

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
        config["kvf"]
    )

    # Save preprocessed data
    # Note: fluxes (dim: time) are at the edges of the timesteps,
    # while states (dim: time2) are at the midpoints and include next midnight
    # so the first state from day 2 will overlap with the last flux from day 1
    filename = f"{date.strftime('%Y-%m-%d')}_fluxes_storages.nc"
    output_path = os.path.join(config["preprocessed_data_folder"], filename)
    xr.Dataset(
        {  # TODO: would be nice to add coordinates and units as well
            "fx_upper": (["time_fluxes", "lat", "lon"], fx_upper, {"units": "m**3"}),
            "fy_upper": (["time_fluxes", "lat", "lon"], fy_upper, {"units": "m**3"}),
            "fx_lower": (["time_fluxes", "lat", "lon"], fx_lower, {"units": "m**3"}),
            "fy_lower": (["time_fluxes", "lat", "lon"], fy_lower, {"units": "m**3"}),
            "s_upper": (["time_states", "lat", "lon"], s_upper, {"units": "m**3"}),
            "s_lower": (["time_states", "lat", "lon"], s_lower, {"units": "m**3"}),
            "f_vert": (["time_fluxes", "lat", "lon"], f_vert, {"units": "m**3"}),
            "evap": (["time_fluxes", "lat", "lon"], evap, {"units": "m**3"}),
            "precip": (["time_fluxes", "lat", "lon"], precip, {"units": "m**3"}),
        },
        coords={
            'time_fluxes': newtime_fluxes,
            'time_states': newtime_states,
            'lat': u.latitude.values,
            'lon': u.longitude.values
        }
    ).to_netcdf(output_path)
