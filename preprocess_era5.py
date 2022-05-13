import calendar
import datetime as dt
import glob
import os
from timeit import default_timer as timer

# Import libraries
import numpy as np
import scipy.io as sio
import xarray as xr
import yaml

from getconstants_pressure_ERA5 import getconstants_pressure_ERA5

from preprocessing import getrefined_new, get_stablefluxes, getFa_Vert

# Read case configuration
with open("cases/era5.yaml") as f:
    config = yaml.safe_load(f)

# Parse input from config file
# Reassignment not strictly needed but improves readability for often used vars
input_folder = config["input_folder"]
name_of_run = config["name_of_run"]
divt = config["divt"]
count_time = config["count_time"]
latnrs = np.arange(config["latnrs"])
lonnrs = np.arange(config["lonnrs"])


# to create datelist
def get_times_daily(startdate, enddate):
    """generate a dictionary with date/times"""
    numdays = enddate - startdate
    dateList = []
    for x in range(0, numdays.days + 1):
        dateList.append(startdate + dt.timedelta(days=x))
    return dateList


months = np.arange(config["start_month"], config["end_month"])
months_length_leap = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
months_length_nonleap = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
years = np.arange(config["start_year"], config["end_year"])

# create datelist
if int(calendar.isleap(years[-1])) == 0:  # no leap year
    datelist = get_times_daily(
        dt.date(years[0], months[0], 1),
        dt.date(years[-1], months[-1], months_length_nonleap[months[-1] - 1]),
    )
else:  # leap
    datelist = get_times_daily(
        dt.date(years[0], months[0], 1),
        dt.date(years[-1], months[-1], months_length_leap[months[-1] - 1]),
    )


def get_3d_data(variable, year, month):
    filename = f"{variable}_ERA5_{year}_{month:02d}_*NH.nc"
    filepath = os.path.join(config["input_folder"], filename)
    files = []
    for file in glob.glob(filepath):
        ds = xr.open_dataset(file, chunks='auto')
        files.append(ds)

    return xr.concat(files, dim="level").sortby("level")[variable]


def get_2d_data(variable, year, month):
    filename = f"{variable}_ERA5_{year}_{month:02d}_NH.nc"
    filepath = os.path.join(config["input_folder"], filename)
    return xr.open_dataset(filepath, chunks='auto')


def get_output_data(year, month, a):
    """Get data for output file."""
    filename = f"{year}-{month:02d}-{a:02d}fluxes_storages.mat"
    save_path = os.path.join(config["interdata_folder"], filename)
    return save_path


def getWandFluxes(date, density_water, g, A_gridcell):
    """Determine the fluxes and states."""

    # Load data
    u = get_3d_data("u", date.year, date.month).sel(time=date.strftime('%Y%m%d'))
    v = get_3d_data("v", date.year, date.month).sel(time=date.strftime('%Y%m%d'))
    q = get_3d_data("q", date.year, date.month).sel(time=date.strftime('%Y%m%d'))
    sp = get_2d_data("sp", date.year, date.month).sel(time=date.strftime('%Y%m%d'))
    levels = q.level

    # Calculate fluxes
    uq = u * q
    vq = v * q

    # Integrate fluxes over original ERA5 layers
    midpoints = 0.5 * (levels.values[1:] + levels.values[:-1])
    uq_midpoints = uq.interp(level=midpoints)
    vq_midpoints = vq.interp(level=midpoints)
    q_midpoints = q.interp(level=midpoints)

    p = levels.broadcast_like(u)  # hPa
    dp_midpoints = p.diff(dim="level")
    dp_midpoints["level"] = midpoints

    # eastward and northward fluxes
    eastward_flux = uq_midpoints * dp_midpoints / g
    northward_flux = vq_midpoints * dp_midpoints / g

    # compute the column water vapor = specific humidity * pressure levels length
    cwv = q_midpoints * dp_midpoints / g

    # total column water vapor, cwv is summed over the vertical [kg/m2]
    tcwv = cwv.sum(dim="level")

    # apply mask so only level in upper or lower layer are masked, based on surface pressure and pressur on boundary
    ##Imme: location of boundary is hereby hard defined at model level 47 which corresponds with about
    P_boundary = 0.72878581 * sp + 7438.803223
    lower_layer = (p.level < sp.sp / 100) & (p.level > P_boundary.sp / 100)
    upper_layer = p.level < P_boundary.sp / 100

    # integrate the fluxes over the upper and lower layer based on the mask
    eastward_flux_lower = eastward_flux.where(lower_layer).sum(dim="level")
    northward_flux_lower = northward_flux.where(lower_layer).sum(dim="level")

    eastward_flux_upper = eastward_flux.where(upper_layer).sum(dim="level")
    northward_flux_upper = northward_flux.where(upper_layer).sum(dim="level")

    W_upper = cwv.where(upper_layer).sum(dim="level")
    W_lower = cwv.where(lower_layer).sum(dim="level")

    vapor_total = W_upper + W_lower

    test0 = (tcwv - vapor_total).sum().values
    print(f"Check calculation water vapor, this value should be zero: {test0}")

    # convert units to water volumes m3
    W_top = W_upper * A_gridcell / density_water  # m3
    W_down = W_lower * A_gridcell / density_water  # m3

    return (
        # TODO: don't load everything in memory here
        cwv.values,
        W_top.values,
        W_down.values,
        eastward_flux_upper.values,
        northward_flux_upper.values,
        eastward_flux_lower.values,
        northward_flux_lower.values,
    )


# Code
def getEP(date, A_gridcell):
    """Load precipitation and evaporation data with units of m3."""
    evap = get_2d_data("e", date.year, date.month)
    precip = get_2d_data("tp", date.year, date.month)

    # calculate volumes
    evap *= A_gridcell  # m3
    precip *= A_gridcell  # m3

    # TODO: don't load everything in memory here
    return evap.values, precip.values

#%% Runtime & Results

start1 = timer()

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
) = getconstants_pressure_ERA5(latnrs, lonnrs, config["land_sea_mask"])

for date in datelist[:2]:
    # start = timer()
    print(("0 = " + str(timer())))
    #    #1 integrate specific humidity to get the (total) column water (vapor) and calculate horizontal moisture fluxes
    cwv, W_top, W_down, Fa_E_top, Fa_N_top, Fa_E_down, Fa_N_down = getWandFluxes(
        date,
        density_water,
        g,
        A_gridcell,
    )
    print(("1,2,3 = " + str(timer())))

    # 4 evaporation and precipitation
    E, P = getEP(date, A_gridcell)
    print(("4 = " + str(timer())))

    # put data on a smaller time step
    (
        Fa_E_top_1,
        Fa_N_top_1,
        Fa_E_down_1,
        Fa_N_down_1,
        E,
        P,
        W_top,
        W_down,
    ) = getrefined_new(
        Fa_E_top,
        Fa_N_top,
        Fa_E_down,
        Fa_N_down,
        W_top,
        W_down,
        E,
        P,
        divt,
        count_time,
        latitude,
        longitude,
    )
    print(("5 = " + str(timer())))

    # change units to m3
    Fa_E_top_m3, Fa_E_down_m3, Fa_N_top_m3, Fa_N_down_m3 = change_units(
        Fa_E_top_1,
        Fa_E_down_1,
        Fa_N_top_1,
        Fa_N_down_1,
        timestep,
        divt,
        L_EW_gridcell,
        density_water,
        L_N_gridcell,
        L_S_gridcell,
        latitude,
    )
    print(("6a = " + str(timer())))

    # stabilize horizontal fluxes
    Fa_E_top, Fa_E_down, Fa_N_top, Fa_N_down = get_stablefluxes(
        Fa_E_top_m3,
        Fa_E_down_m3,
        Fa_N_top_m3,
        Fa_N_down_m3,
        timestep,
        divt,
        L_EW_gridcell,
        density_water,
        L_N_gridcell,
        L_S_gridcell,
        latitude,
    )
    print(("6b = " + str(timer())))

    # determine the vertical moisture flux
    Fa_Vert_raw, Fa_Vert, residual_down, residual_top = getFa_Vert(
        Fa_E_top,
        Fa_E_down,
        Fa_N_top,
        Fa_N_down,
        E,
        P,
        W_top,
        W_down,
        divt,
        count_time,
        latitude,
        longitude,
    )
    print(("7 = " + str(timer())))

    # np.savez_compressed(get_output_path(year, month, a), E=E, P=P, Fa_E_top=Fa_E_top, Fa_N_top= Fa_N_top, Fa_E_down=Fa_E_down, Fa_N_down=Fa_N_down, W_down=W_down, W_top=W_top, residual_top=residual_top, residual_down=residual_down, Fa_Vert=Fa_Vert) # save as .npy file
    sio.savemat(
        get_output_path(year, month, a),
        {
            "Fa_E_top": Fa_E_top,
            "Fa_N_top": Fa_N_top,
            "Fa_E_down": Fa_E_down,
            "Fa_N_down": Fa_N_down,
            "E": E,
            "P": P,
            "W_top": W_top,
            "W_down": W_down,
            "Fa_Vert": Fa_Vert,
        },
        do_compression=True,
    )  # save as mat file

    end = timer()
    print(
        "Runtime fluxes_and_storages for day "
        + str(a)
        + " in year "
        + str(yearnumber)
        + " is",
        (end - start),
        " seconds.",
    )
end1 = timer()
print("The total runtime is", (end1 - start1), " seconds.")
