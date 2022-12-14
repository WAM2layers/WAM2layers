import datetime as dt
import os

import numpy as np
import pandas as pd
import xarray as xr
import yaml

from .shared import (get_grid_info, get_new_target_levels, get_stable_fluxes,
                     get_vertical_transport, interpolate, join_levels,
                     repeat_upper_level, resample)

# Set constants
g = 9.80665  # [m/s2]
density_water = 1000  # [kg/m3]

# Read case configuration
with open("cases/ec-earth.yaml") as f:
    config = yaml.safe_load(f)

# Parse input from config file
# Reassignment not strictly needed but improves readability for often used vars
input_folder = config["input_folder"]
name_of_run = config["name_of_run"]
divt = config["divt"]
count_time = config["count_time"]
latnrs = np.arange(config["latnrs"])
lonnrs = np.arange(config["lonnrs"])


def _get_input_data(variable, date, latnrs, lonnrs):
    """Get input data for variable."""
    filename = f"{name_of_run}{variable}_{date.year}{date.month:02d}_NH.nc"
    filepath = os.path.join(config["input_folder"], filename)
    return (
        xr.open_dataset(filepath)
        .sel(time=date.strftime("%Y%m%d"))
        .isel(lat=latnrs, lon=lonnrs)
    )


def get_input_data(variable, date, latnrs, lonnrs):
    """Workaround for keeping original timestamps crossing midnight."""
    day1 = _get_input_data(variable, date, latnrs, lonnrs)
    day2 = _get_input_data(variable, date + dt.timedelta(days=1), latnrs, lonnrs)
    days = xr.concat([day1, day2], dim="time")
    if len(days.time) < 9:
        return days.isel(time=slice(0, 5)).isel(lat=latnrs, lon=lonnrs)
    else:
        return days.isel(time=slice(0, 10)).isel(lat=latnrs, lon=lonnrs)


datelist = pd.date_range(
    start=config["start_date"], end=config["end_date"], freq="d", inclusive="left"
)

for date in datelist:
    print(date)

    # Load data
    time_3h = slice(1, None, 2)
    times_lnsp = slice(None, None, 2)  # TODO should be slice(1, None, 2) ??

    q = get_input_data("Q", date, latnrs, lonnrs).Q
    u = get_input_data("U", date, latnrs, lonnrs).U
    v = get_input_data("V", date, latnrs, lonnrs).V
    q2m = get_input_data("Q2M", date, latnrs, lonnrs).Q2M.isel(time=time_3h)
    lnsp = get_input_data("LNSP", date, latnrs, lonnrs).LNSP.isel(time=times_lnsp)
    u10 = get_input_data("U10", date, latnrs, lonnrs).U10M.isel(time=time_3h)
    v10 = get_input_data("V10", date, latnrs, lonnrs).V10M.isel(time=time_3h)

    sp = np.exp(lnsp)  # log(sp) --> sp (Pa)

    # Create pressure array with the same dimensions as u, q, and v
    p = u.lev.broadcast_like(u)

    # Get the most out of the data we have
    u = join_levels(u, u10)  # bottom level will be set to 1000 00 Pa
    v = join_levels(v, v10)
    q = join_levels(q, q2m)
    p = join_levels(p, sp.squeeze())

    u = repeat_upper_level(u)  # top level will be set to 100 00 Pa
    v = repeat_upper_level(v)
    q = repeat_upper_level(q, fill_value=0)
    p = repeat_upper_level(p, fill_value=0)

    # Mask data where the pressure level is higher than sp - 1000
    # as there is no valid data at those points
    # and we don't need a second layer that is very close to the surface
    mask = p > sp.values - 1000
    mask[:, 0, :, :] = False  # don't mask bottom (surface pressure values)
    mask[:, -1, :, :] = False  # don't mask top ("ghost cells"?)

    # TODO convert to nan instead of numpy masked array?
    # u_masked = u.where(mask)
    u = np.ma.masked_array(u, mask=mask)
    v = np.ma.masked_array(v, mask=mask)
    q = np.ma.masked_array(q, mask=mask)
    p = np.ma.masked_array(p, mask=mask)

    ####
    # Get grid info
    dummy = xr.open_dataset(config["land_sea_mask"])
    lat = dummy.LAT.values[444:711][::-1]  # [degrees north]
    lon = dummy.XAS.values[934:1378][::-1]  # [degrees east]
    a_gridcell, l_ew_gridcell, l_mid_gridcell = get_grid_info(lat, lon)

    # Reverse engineered ~ model level 47 which corresponds with about ~800 hPa
    p_boundary = 0.72878581 * sp.values + 7438.803223
    new_pressure_levels = get_new_target_levels(sp, p_boundary, n_levels=40)

    print("before interpolation loop", dt.datetime.now().time())
    uq_boundaries = interpolate(u * q, p, new_pressure_levels, type="cubic")
    vq_boundaries = interpolate(v * q, p, new_pressure_levels, type="cubic")
    q_boundaries = interpolate(q, p, new_pressure_levels, type="linear")

    # Integrate specific humidity to get the (total) column water (vapor) and calculate horizontal moisture fluxes
    q_midpoints = 0.5 * (q_boundaries[:, 1:, :, :] + q_boundaries[:, :-1, :, :])
    uq_midpoints = 0.5 * (uq_boundaries[:, 1:, :, :] + uq_boundaries[:, :-1, :, :])
    vq_midpoints = 0.5 * (vq_boundaries[:, 1:, :, :] + vq_boundaries[:, :-1, :, :])
    # for p we do not calculate the mean pressure but the pressure difference between two levels!
    p_diff = np.maximum(
        0, new_pressure_levels[:, :-1, :, :] - new_pressure_levels[:, 1:, :, :]
    )  # the maximum statement is necessary to avoid negative humidity values

    # eastward and northward fluxes
    fa_e_p = uq_midpoints * p_diff / g
    fa_n_p = vq_midpoints * p_diff / g

    column_water_vapor = q_midpoints * p_diff / g  # [kg/m2]

    # summed over the vertical
    total_column_water_vapor = np.squeeze(np.sum(column_water_vapor, 1))

    mask = np.where(new_pressure_levels > p_boundary, 1.0, 0.0)

    vapor_lower = np.sum(mask[:, :-1, :, :] * q_midpoints * p_diff / g, axis=1)
    vapor_upper = np.sum((1 - mask[:, :-1, :, :]) * q_midpoints * p_diff / g, axis=1)

    fa_e_lower = np.sum(mask[:, :-1, :, :] * fa_e_p, axis=1)  # kg*m-1*s-1
    fa_n_lower = np.sum(mask[:, :-1, :, :] * fa_n_p, axis=1)  # kg*m-1*s-1
    fa_e_upper = np.sum((1 - mask[:, :-1, :, :]) * fa_e_p, axis=1)  # kg*m-1*s-1
    fa_n_upper = np.sum((1 - mask[:, :-1, :, :]) * fa_n_p, axis=1)  # kg*m-1*s-1

    vapor_total = vapor_upper + vapor_lower

    check = np.sum(total_column_water_vapor - vapor_total)
    print(f"Check calculation water vapor, this value should be zero: {check}")

    # water volumes
    w_upper = vapor_upper * a_gridcell[None, ...] / density_water  # m3
    w_lower = vapor_lower * a_gridcell[None, ...] / density_water  # m3

    ## Get E and P
    # (accumulated after the forecast at 00.00 and 12.00 by steps of 3 hours in time
    evaporation = get_input_data("EVAP", date, latnrs, lonnrs).E  # m
    precipitation = get_input_data("TP", date, latnrs, lonnrs).TP  # m

    # change sign convention to all positive, transfer negative (originally positive) values of evap to precip
    precipitation = np.maximum(precipitation + np.maximum(evaporation, 0), 0)
    evaporation = np.abs(np.minimum(evaporation, 0))

    # calculate volumes
    evap = (evaporation * a_gridcell).values
    precip = (precipitation * a_gridcell).values

    # put data on a smaller time step
    evap = resample(evap, divt / 2, count_time, method="bfill")[:-1]
    precip = resample(precip, divt / 2, count_time, method="bfill")[:-1]
    w_upper = resample(w_upper, divt, count_time, method="interp")
    w_lower = resample(w_lower, divt, count_time, method="interp")
    fa_e_upper = resample(fa_e_upper, divt, count_time, method="interp")[:-1]
    fa_e_lower = resample(fa_e_lower, divt, count_time, method="interp")[:-1]
    fa_n_upper = resample(fa_n_upper, divt, count_time, method="interp")[:-1]
    fa_n_lower = resample(fa_n_lower, divt, count_time, method="interp")[:-1]

    # convert to m3   [kg*m^-1 * s^-1 * s * m * kg^-1 * m^3] = [m3]
    total_seconds = config["timestep"] / config["divt"]
    fa_e_upper *= total_seconds * (l_ew_gridcell / density_water)
    fa_e_lower *= total_seconds * (l_ew_gridcell / density_water)
    fa_n_upper *= total_seconds * (l_mid_gridcell[None, :, None] / density_water)
    fa_n_lower *= total_seconds * (l_mid_gridcell[None, :, None] / density_water)

    # stabilize horizontal fluxes
    fa_e_upper, fa_n_upper = get_stable_fluxes(fa_e_upper, fa_n_upper, w_upper)
    fa_e_lower, fa_n_lower = get_stable_fluxes(fa_e_lower, fa_n_lower, w_lower)

    # determine the vertical moisture flux
    fa_vert = get_vertical_transport(
        fa_e_upper, fa_e_lower, fa_n_upper, fa_n_lower, evap, precip, w_upper, w_lower
    )

    # Save preprocessed data
    filename = f"{date.strftime('%Y-%m-%d')}_fluxes_storages.nc"
    output_path = os.path.join(config["interdata_folder"], filename)
    xr.Dataset(
        {  # TODO: would be nice to add coordinates and units as well
            "fa_e_upper": (["time", "lat", "lon"], fa_e_upper),
            "fa_n_upper": (["time", "lat", "lon"], fa_n_upper),
            "fa_e_lower": (["time", "lat", "lon"], fa_e_lower),
            "fa_n_lower": (["time", "lat", "lon"], fa_n_lower),
            "evap": (["time", "lat", "lon"], evap),
            "precip": (["time", "lat", "lon"], precip),
            "w_upper": (["time2", "lat", "lon"], w_upper),  # note different time
            "w_lower": (["time2", "lat", "lon"], w_lower),  # note different time
            "fa_vert": (["time", "lat", "lon"], fa_vert),
        }
    ).to_netcdf(output_path)
