import os

# import numpy as np
import pandas as pd
import xarray as xr
import yaml

from getconstants_pressure_ERA5 import getconstants_pressure_ERA5
from preprocessing import get_stablefluxes, getFa_Vert

# Read case configuration
with open("cases/era5_2013.yaml") as f:
    config = yaml.safe_load(f)

# Parse input from config file
# Reassignment not strictly needed but improves readability for often used vars
input_folder = config["input_folder"]
name_of_run = config["name_of_run"]
divt = config["divt"]
count_time = config["count_time"]

datelist = pd.date_range(
    start=config["start_date"], end=config["end_date"], freq="d", inclusive="left"
)


def get_data(variable, year, month):
    # filename = f"{variable}_ERA5_{year}_{month:02d}_NH.nc"
    filename = f"FloodCase_{year}{month:02d}_{variable}.nc"
    filepath = os.path.join(config["input_folder"], filename)
    return xr.open_dataset(filepath)[variable]


def get_output_path(year, month, day, extension=".nc"):
    """Get data for output file."""
    filename = f"{year}-{month:02d}-{day:02d}fluxes_storages{extension}"
    save_path = os.path.join(config["interdata_folder"], filename)
    return save_path

# obtain the constants
(
    latitude,
    longitude,
    lsm,
    g,
    density_water,
    timestep,
    A_gridcell,
    L_N_gridcell,
    L_S_gridcell,
    L_EW_gridcell,
    gridcell,
) = getconstants_pressure_ERA5(config["land_sea_mask"])

for date in datelist:
    print(date)

    year = date.year
    month = date.month
    day = date.day

    # Load data
    u = get_data("u", year, month).sel(time=date.strftime('%Y%m%d'))
    v = get_data("v", year, month).sel(time=date.strftime('%Y%m%d'))
    q = get_data("q", year, month).sel(time=date.strftime('%Y%m%d'))
    sp = get_data("sp", year, month).sel(time=date.strftime('%Y%m%d'))
    evap = get_data("e", year, month).sel(time=date.strftime('%Y%m%d'))
    cp = get_data("cp", year, month).sel(time=date.strftime('%Y%m%d'))
    lsp = get_data("lsp", year, month).sel(time=date.strftime('%Y%m%d'))
    precip = cp + lsp

    # Calculate volumes
    evap *= A_gridcell  # m3
    precip *= A_gridcell  # m3

    # Create pressure array
    levels = q.level
    p = levels.broadcast_like(u)  # hPa

    # Interpolate to new levels
    edges = 0.5 * (levels.values[1:] + levels.values[:-1])
    u = u.interp(level=edges)
    v = v.interp(level=edges)
    q = q.interp(level=edges)

    # Calculate pressure jump
    dp = p.diff(dim="level")
    dp["level"] = edges

    # Determine the fluxes and states
    uq = u*q * dp / g  # eastward moisture flux
    vq = v*q * dp / g  # northward moisture flux
    cwv = q * dp / g * A_gridcell / density_water  # column water vapor (m3)

    # Split in 2 layers
    P_boundary = 0.72878581 * sp + 7438.803223
    lower_layer = (dp.level < sp / 100) & (dp.level > P_boundary / 100)
    upper_layer = dp.level < P_boundary / 100

    # Integrate fluxes and state
    uq_lower = uq.where(lower_layer).sum(dim="level")
    vq_lower = vq.where(lower_layer).sum(dim="level")
    w_lower = cwv.where(lower_layer).sum(dim="level")

    uq_upper = uq.where(upper_layer).sum(dim="level")
    vq_upper = vq.where(upper_layer).sum(dim="level")
    w_upper = cwv.where(upper_layer).sum(dim="level")

    print(
        "Check calculation water vapor, this value should be zero:",
        (cwv.sum(dim="level") - (w_upper + w_lower)).sum().values
    )

    # Change units to m3
    density_water = 1000  # [kg/m3]
    L_mid_gridcell = 0.5 * (L_N_gridcell + L_S_gridcell)
    uq_upper *= (timestep / divt) * (L_EW_gridcell / density_water)
    uq_lower *= (timestep / divt) * (L_EW_gridcell / density_water)
    vq_upper *= (timestep / divt) * (L_mid_gridcell[None, :, None] / density_water)
    vq_lower *= (timestep / divt) * (L_mid_gridcell[None, :, None] / density_water)

    # Put data on a smaller time step
    time = w_upper.time.values
    newtime = pd.date_range(time[0], time[-1], freq='15Min')
    w_upper = w_upper.interp(time=newtime).values
    w_lower = w_lower.interp(time=newtime).values

    # Put fluxes on the edges instead of midpoints
    newtime = newtime[:-1] + pd.Timedelta('15Min') / 2
    uq_upper = uq_upper.interp(time=newtime).values
    vq_upper = vq_upper.interp(time=newtime).values

    uq_lower = uq_lower.interp(time=newtime).values
    vq_lower = vq_lower.interp(time=newtime).values

    precip = (precip.reindex(time=newtime, method='bfill') / 4).values
    evap = (evap.reindex(time=newtime, method='bfill') / 4).values

    # Stabilize horizontal fluxes
    uq_upper, uq_upper = get_stablefluxes(uq_upper, vq_upper, w_upper)
    uq_lower, uq_lower = get_stablefluxes(uq_lower, vq_lower, w_lower)

    # Determine the vertical moisture flux
    Fa_Vert = getFa_Vert(
        uq_upper,
        uq_lower,
        vq_upper,
        vq_lower,
        evap,
        precip,
        w_upper,
        w_lower,
        )

    # Save preprocessed data
    xr.Dataset(
        {  # TODO: would be nice to add coordinates and units as well
            "uq_upper": (["time", "lat", "lon"], uq_upper),
            "vq_upper": (["time", "lat", "lon"], vq_upper),
            "uq_lower": (["time", "lat", "lon"], uq_lower),
            "vq_lower": (["time", "lat", "lon"], vq_lower),
            "w_upper": (["time2", "lat", "lon"], w_upper),  # note different time
            "w_lower": (["time2", "lat", "lon"], w_lower),  # note different time
            "Fa_Vert": (["time", "lat", "lon"], Fa_Vert),
            "E": (["time", "lat", "lon"], evap),
            "P": (["time", "lat", "lon"], precip),
        }
    ).to_netcdf(get_output_path(year, month, day, extension=".nc"))
