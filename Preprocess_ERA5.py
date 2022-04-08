# This script is similar as the Fluxes_and_States_Masterscript from WAM-2layers from Ruud van der Ent, except that it threads atmospheric data at five pressure levels instead of model levels

# Includes a spline vertical interpolation for data on pressure levels from the EC-Earth climate model (Hazeleger et al., 2010)
# In this case the atmospheric data is provided on the following five pressure levels: 850 hPa, 700 hPa, 500 hPa, 300 hPa, 200 hPa
# Includes a linear interpolation of the moisture fluxes over time (in def getrefined_new)
# We have implemented a datelist function so the model can run for multiple years without having problems with leap years

# My input data are monthly files with hourly (surface variables) and hourly data (atmospheric variables)

import calendar
import datetime as dt
import os
import sys
from datetime import timedelta
from timeit import default_timer as timer

import matplotlib.pyplot as plt
# Import libraries
import numpy as np
import scipy.io as sio
import timer
import xarray as xr
import yaml
from netCDF4 import Dataset
from scipy import interpolate
from scipy.interpolate import interp1d

from getconstants_pressure_ERA5 import getconstants_pressure_ERA5

# Read case configuration
with open("cases/example.yaml") as f:
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


def get_input_data(variable, year, month):
    """Get input data for variable."""
    filename = f"{name_of_run}{variable}_ERA5_{year}_{month:02d}_trop_NH.nc"
    filepath = os.path.join(config["input_folder"], filename)
    return Dataset(filepath, "r")


def get_3d_data(variable, year, month):
    filename = f"{variable}_ERA5_{year}_{month:02d}_*NH.nc"
    filepath = os.path.join(config["input_folder"], filename)
    files = []
    for file in glob.glob(filepath):
        ds = xr.open_dataset(file)
        files.append(ds)

    return xr.concat(files, dim='level').sortby('level')

def get_2d_data(variable, year, month):
    filename = f"{variable}_ERA5_{year}_{month:02d}_*NH.nc"
    filepath = os.path.join(config["input_folder"], filename)
    return xr.open_mfdataset(filepath)


def get_input_data_next_month(variable, year, month):
    """Get input data for next month(?)"""
    if month == 12:
        return get_input_data(variable, year + 1, 1)
    else:
        return get_input_data(variable, year, month + 1)


def get_output_data(year, month, a):
    """Get data for output file."""
    filename = f"{year}-{month:02d}-{a:02d}fluxes_storages.mat"
    save_path = os.path.join(config["interdata_folder"], filename)
    return save_path


# Code (no need to look at this for running)

# Determine the fluxes and states
# In this defintion the vertical spline interpolation is performed to determine the moisture fluxes for the two layers at each grid cell
def getWandFluxes(
    latnrs,
    lonnrs,
    date,
    density_water,
    latitude,
    longitude,
    g,
    A_gridcell,
):

    q = get_3d_data("q", date.year, date.month).sel(time=date.strftime('%Y%m%d'))
    time = q.time
    levels = q.level

    u = get_3d_data("u", date.year, date.month).sel(time=date.strftime('%Y%m%d'))
    v = get_3d_data("v", date.year, date.month).sel(time=date.strftime('%Y%m%d'))
    sp = get_2d_data("sp", date.year, date.month).sel(time=date.strftime('%Y%m%d'))

    print("Data is loaded", dt.datetime.now().time())
    #time = [0, 1, 2, 3, 4] # do not hardcode this
    # If you do not have 5 pressure levels but much more, as in our case, then calculate 
    # intervals directly from the real data
    #but then also take surface pressure into account and remove that

    ##Imme: location of boundary is hereby hard defined at model level 47 which corresponds with about
    P_boundary = 0.72878581 * sp + 7438.803223
    #dp = (sp - 20000.0) / (intervals_regular - 1) # this has to change, here you determine the difference between the pressur levels


    uq = u*q
    vq = v*q

    q_midpoints = 0.5 * (q.isel(level=slice(1, None)) + q.isel(level=slice(None, -1)))
    uq_midpoints = 0.5 * (uq.isel(level=slice(1, None)) + uq.isel(level=slice(None, -1)))
    vq_midpoints = 0.5 * (vq.isel(level=slice(1, None)) + vq.isel(level=slice(None, -1)))

    midpoints = 0.5 * (levels.values[1:] + levels.values[:-1])
    

    u_lower = u.where((u.level < sp.sp/100) & (u.level > P_boundary.sp/100))
    v_lower = v.where((v.level < sp.sp/100) & (v.level > P_boundary.sp/100))
    q_lower = q.where((q.level < sp.sp/100) & (q.level > P_boundary.sp/100))
    
    u_upper = u.where((u.level < P_boundary.sp/100))
    v_upper = v.where((v.level < P_boundary.sp/100))
    q_upper = q.where((q.level < P_boundary.sp/100))

    # p.diff(dim='level')
    # levels.broadcast_like(u)
    # p = sp.broadcast_like(u)

    # eerst alles interpoleren naar de tussen levels
    # ook p interpoleren naar tussen levels
    # dp uitrekenen
    # dan u en q en zo vermenidvuldigen en masken
    levels.diff(dim='level')

    dp = 
    #TODO: how to deal with lowest layers, interpolate to take into account level below surface to have more information
    p_maxmin = np.zeros((time, intervals_regular + 2, len(latitude), len(longitude)))
    p_maxmin[:, 1:-1, :, :] = (
        sp[:, np.newaxis, :, :]
        - dp[:, np.newaxis, :, :]
        * np.arange(0, intervals_regular)[np.newaxis, :, np.newaxis, np.newaxis]
    )

    mask = np.where(p_maxmin > P_boundary[:, np.newaxis, :, :], 1.0, 0.0)
    mask[:, 0, :, :] = 1.0  # bottom value is always 1

    p_maxmin[:, :-1, :, :] = (
        mask[:, 1:, :, :] * p_maxmin[:, 1:, :, :]
        + (1 - mask[:, 1:, :, :]) * p_maxmin[:, :-1, :, :]
    )
    p_maxmin[:, 1:, :, :] = np.where(
        p_maxmin[:, :-1, :, :] == p_maxmin[:, 1:, :, :],
        P_boundary[:, np.newaxis, :, :],
        p_maxmin[:, 1:, :, :],
    )

    del (dp, mask)
    print(
        "after p_maxmin and now add surface and atmosphere together for u,q,v",
        dt.datetime.now().time(),
    )

    levelist = np.squeeze(get_input_data("U", year, month).variables["lev"])  # Pa
    p = np.zeros((time, levelist.size + 2, len(latitude), len(longitude)))
    p[:, 1:-1, :, :] = levelist[np.newaxis, :, np.newaxis, np.newaxis]
    p[:, 0, :, :] = sp

    u_total = np.zeros((time, levelist.size + 2, len(latitude), len(longitude)))
    u_total[:, 1:-1, :, :] = u
    u_total[:, 0, :, :] = u10 # do not add the 10m with ERA5 data that is not there
    u_total[:, -1, :, :] = u[:, -1, :, :]

    v_total = np.zeros((time, levelist.size + 2, len(latitude), len(longitude)))
    v_total[:, 1:-1, :, :] = v
    v_total[:, 0, :, :] = v10 # do not add the 10m with ERA5 data that is not there
    v_total[:, -1, :, :] = v[:, -1, :, :]

    q_total = np.zeros((time, levelist.size + 2, len(latitude), len(longitude)))
    q_total[:, 1:-1, :, :] = q
    q_total[:, 0, :, :] = q2m # do not add the 2m q with ERA5 data that is not there

    mask = np.ones(u_total.shape, dtype=np.bool)
    mask[:, 1:-1, :, :] = levelist[np.newaxis, :, np.newaxis, np.newaxis] < (
        sp[:, np.newaxis, :, :] - 1000.0
    )  # Pa

    u_masked = np.ma.masked_array(u_total, mask=~mask)
    v_masked = np.ma.masked_array(v_total, mask=~mask)
    q_masked = np.ma.masked_array(q_total, mask=~mask)
    p_masked = np.ma.masked_array(p, mask=~mask)

    del (u_total, v_total, q_total, p, u, v, q, u10, v10, q2m, sp)

    print("before interpolation loop", dt.datetime.now().time())

    uq_maxmin = np.zeros((time, intervals_regular + 2, len(latitude), len(longitude)))
    vq_maxmin = np.zeros((time, intervals_regular + 2, len(latitude), len(longitude)))
    q_maxmin = np.zeros((time, intervals_regular + 2, len(latitude), len(longitude)))

    #only do this when there are not enought pressure levels available in the dataset
    # otherwiste just a linear interpolation.
    for t in range(time):  # loop over timesteps
        for i in range(len(latitude)):  # loop over latitude
            for j in range(len(longitude)):  # loop over longitude
                pp = p_masked[t, :, i, j]
                uu = u_masked[t, :, i, j]
                vv = v_masked[t, :, i, j]
                qq = q_masked[t, :, i, j]

                pp = pp[~pp.mask]
                uu = uu[~uu.mask]
                vv = vv[~vv.mask]
                qq = qq[~qq.mask]

                f_uq = interp1d(pp, uu * qq, "cubic")  # spline interpolation
                uq_maxmin[t, :, i, j] = f_uq(
                    p_maxmin[t, :, i, j]
                )  # spline interpolation

                f_vq = interp1d(pp, vv * qq, "cubic")  # spline interpolation
                vq_maxmin[t, :, i, j] = f_vq(
                    p_maxmin[t, :, i, j]
                )  # spline interpolation

                f_q = interp1d(pp, qq)  # linear interpolation
                q_maxmin[t, :, i, j] = f_q(p_maxmin[t, :, i, j])  # linear interpolation

    del (u_masked, v_masked, q_masked, p_masked, mask, f_uq, f_vq, f_q)
    print("after interpolation loop", dt.datetime.now().time())

    # pressure between full levels
    P_between = np.maximum(
        0, p_maxmin[:, :-1, :, :] - p_maxmin[:, 1:, :, :]
    )  # the maximum statement is necessary to avoid negative humidity values
    # Imme: in P_between you do not calculate the pressure between two levels but the pressure difference between two levels!!!
    q_between = 0.5 * (q_maxmin[:, 1:, :, :] + q_maxmin[:, :-1, :, :])
    uq_between = 0.5 * (uq_maxmin[:, 1:, :, :] + uq_maxmin[:, :-1, :, :])
    vq_between = 0.5 * (vq_maxmin[:, 1:, :, :] + vq_maxmin[:, :-1, :, :])

    # eastward and northward fluxes
    Fa_E_p = uq_between * P_between / g
    Fa_N_p = vq_between * P_between / g

    # compute the column water vapor
    cwv = (
        q_between * P_between / g
    )  # column water vapor = specific humidity * pressure levels length / g [kg/m2]
    # make tcwv vector
    tcwv = np.squeeze(
        np.sum(cwv, 1)
    )  # total column water vapor, cwv is summed over the vertical [kg/m2]

    # use mask
    mask = np.where(p_maxmin > P_boundary[:, np.newaxis, :, :], 1.0, 0.0)

    vapor_down = np.sum(mask[:, :-1, :, :] * q_between * P_between / g, axis=1)
    vapor_top = np.sum((1 - mask[:, :-1, :, :]) * q_between * P_between / g, axis=1)

    Fa_E_down = np.sum(mask[:, :-1, :, :] * Fa_E_p, axis=1)  # kg*m-1*s-1
    Fa_N_down = np.sum(mask[:, :-1, :, :] * Fa_N_p, axis=1)  # kg*m-1*s-1
    Fa_E_top = np.sum((1 - mask[:, :-1, :, :]) * Fa_E_p, axis=1)  # kg*m-1*s-1
    Fa_N_top = np.sum((1 - mask[:, :-1, :, :]) * Fa_N_p, axis=1)  # kg*m-1*s-1

    vapor_total = vapor_top + vapor_down

    # check whether the next calculation results in all zeros
    test0 = tcwv - vapor_total
    print(
        (
            "check calculation water vapor, this value should be zero: "
            + str(np.sum(test0))
        )
    )

    # put A_gridcell on a 3D grid
    A_gridcell2D = np.tile(A_gridcell, [1, len(longitude)])
    A_gridcell_1_2D = np.reshape(A_gridcell2D, [1, len(latitude), len(longitude)])
    A_gridcell_plus3D = np.tile(A_gridcell_1_2D, [count_time + 1, 1, 1])

    # water volumes
    W_top = vapor_top * A_gridcell_plus3D / density_water  # m3
    W_down = vapor_down * A_gridcell_plus3D / density_water  # m3

    return cwv, W_top, W_down, Fa_E_top, Fa_N_top, Fa_E_down, Fa_N_down

# Code
def getEP(
    latnrs, lonnrs, year, month, begin_time, count_time, latitude, longitude, A_gridcell
):
    # hourly data so 24 steps per day

    print(begin_time,count_time)
    evaporation = get_input_data("e", year, month).variables["e"][
        begin_time: (begin_time + count_time ), latnrs, lonnrs
    ]  # m
    precipitation = get_input_data("tp", year, month).variables["tp"][
        begin_time : (begin_time + count_time ), latnrs, lonnrs
    ]  # m

    # delete and transfer negative values, change sign convention to all positive
    precipitation = np.reshape(
        np.maximum(
            np.reshape(precipitation, (np.size(precipitation)))
            + np.maximum(np.reshape(evaporation, (np.size(evaporation))), 0.0),
            0.0,
        ),
        (np.int(count_time), len(latitude), len(longitude)),
    )
    evaporation = np.reshape(
        np.abs(np.minimum(np.reshape(evaporation, (np.size(evaporation))), 0.0)),
        (np.int(count_time), len(latitude), len(longitude)),
    )

    # calculate volumes
    A_gridcell2D = np.tile(A_gridcell, [1, len(longitude)])
    A_gridcell_1_2D = np.reshape(A_gridcell2D, [1, len(latitude), len(longitude)])
    A_gridcell_max3D = np.tile(A_gridcell_1_2D, [count_time, 1, 1])

    E = evaporation * A_gridcell_max3D
    P = precipitation * A_gridcell_max3D

    return E, P


# Code

# within this new definition of refined I do a linear interpolation over time of my fluxes
def getrefined_new(
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
):

    # This definition refines the timestep of the data
    # Imme: change the timesteps from 6-hourly and 3-hourly to 96 timesteps a day

    # for 3 hourly information
    divt2 = divt / 2.0
    oddvector2 = np.zeros((1, np.int(count_time * 2 * divt2)))
    partvector2 = np.zeros((1, np.int(count_time * 2 * divt2)))
    da = np.arange(1, divt2)

    for o in np.arange(0, np.int(count_time * 2 * divt2), np.int(divt2)):
        for i in range(len(da)):
            oddvector2[0, o + i] = (divt2 - da[i]) / divt2
            partvector2[0, o + i + 1] = da[i] / divt2

    E_small = np.nan * np.zeros(
        (np.int(count_time * 2 * divt2), len(latitude), len(longitude))
    )
    for t in range(1, np.int(count_time * 2 * divt2) + 1):
        E_small[t - 1] = (1.0 / divt2) * E[np.int(t / divt2 + oddvector2[0, t - 1] - 1)]
    E = E_small

    P_small = np.nan * np.zeros(
        (np.int(count_time * 2 * divt2), len(latitude), len(longitude))
    )
    for t in range(1, np.int(count_time * 2 * divt2) + 1):
        P_small[t - 1] = (1.0 / divt2) * P[np.int(t / divt2 + oddvector2[0, t - 1] - 1)]
    P = P_small

    # for 6 hourly info
    oddvector = np.zeros((1, np.int(count_time * divt)))
    partvector = np.zeros((1, np.int(count_time * divt)))
    da = np.arange(1, divt)
    divt = np.float(divt)
    for o in np.arange(0, np.int(count_time * divt), np.int(divt)):
        for i in range(len(da)):
            oddvector[0, o + i] = (divt - da[i]) / divt
            partvector[0, o + i + 1] = da[i] / divt

    W_top_small = np.nan * np.zeros(
        (np.int(count_time * divt + 1), len(latitude), len(longitude))
    )
    W_down_small = np.nan * np.zeros(
        (np.int(count_time * divt + 1), len(latitude), len(longitude))
    )

    Fa_E_down_small = np.nan * np.zeros(
        (np.int(count_time * divt), len(latitude), len(longitude))
    )
    Fa_N_down_small = np.nan * np.zeros(
        (np.int(count_time * divt), len(latitude), len(longitude))
    )
    Fa_E_top_small = np.nan * np.zeros(
        (np.int(count_time * divt), len(latitude), len(longitude))
    )
    Fa_N_top_small = np.nan * np.zeros(
        (np.int(count_time * divt), len(latitude), len(longitude))
    )

    for t in range(1, np.int(count_time * divt) + 1):
        W_top_small[t - 1] = W_top[
            np.int(t / divt + oddvector[0, t - 1] - 1)
        ] + partvector[0, t - 1] * (
            W_top[np.int(t / divt + oddvector[0, t - 1])]
            - W_top[np.int(t / divt + oddvector[0, t - 1] - 1)]
        )
        W_top_small[-1] = W_top[-1]
        W_down_small[t - 1] = W_down[
            np.int(t / divt + oddvector[0, t - 1] - 1)
        ] + partvector[0, t - 1] * (
            W_down[np.int(t / divt + oddvector[0, t - 1])]
            - W_down[np.int(t / divt + oddvector[0, t - 1] - 1)]
        )
        W_down_small[-1] = W_down[-1]

        Fa_E_down_small[t - 1] = Fa_E_down[
            np.int(t / divt + oddvector[0, t - 1] - 1)
        ] + partvector[0, t - 1] * (
            Fa_E_down[np.int(t / divt + oddvector[0, t - 1])]
            - Fa_E_down[np.int(t / divt + oddvector[0, t - 1] - 1)]
        )
        Fa_N_down_small[t - 1] = Fa_N_down[
            np.int(t / divt + oddvector[0, t - 1] - 1)
        ] + partvector[0, t - 1] * (
            Fa_N_down[np.int(t / divt + oddvector[0, t - 1])]
            - Fa_N_down[np.int(t / divt + oddvector[0, t - 1] - 1)]
        )
        Fa_E_top_small[t - 1] = Fa_E_top[
            np.int(t / divt + oddvector[0, t - 1] - 1)
        ] + partvector[0, t - 1] * (
            Fa_E_top[np.int(t / divt + oddvector[0, t - 1])]
            - Fa_E_top[np.int(t / divt + oddvector[0, t - 1] - 1)]
        )
        Fa_N_top_small[t - 1] = Fa_N_top[
            np.int(t / divt + oddvector[0, t - 1] - 1)
        ] + partvector[0, t - 1] * (
            Fa_N_top[np.int(t / divt + oddvector[0, t - 1])]
            - Fa_N_top[np.int(t / divt + oddvector[0, t - 1] - 1)]
        )

    W_top = W_top_small
    W_down = W_down_small
    Fa_E_down = Fa_E_down_small
    Fa_N_down = Fa_N_down_small
    Fa_E_top = Fa_E_top_small
    Fa_N_top = Fa_N_top_small

    return Fa_E_top, Fa_N_top, Fa_E_down, Fa_N_down, E, P, W_top, W_down


# Code
def change_units(
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
):

    # redefine according to units
    Fa_E_top_kgpmps = Fa_E_top_1
    Fa_E_down_kgpmps = Fa_E_down_1
    Fa_N_top_kgpmps = Fa_N_top_1
    Fa_N_down_kgpmps = Fa_N_down_1

    # convert to m3
    Fa_E_top_m3 = (
        Fa_E_top_kgpmps * timestep / np.float(divt) * L_EW_gridcell / density_water
    )  # [kg*m^-1*s^-1*s*m*kg^-1*m^3]=[m3]
    Fa_E_down_m3 = (
        Fa_E_down_kgpmps * timestep / np.float(divt) * L_EW_gridcell / density_water
    )  # [s*m*kg*m^-1*s^-1*kg^-1*m^3]=[m3]

    Fa_N_top_swap = np.zeros(
        (len(latitude), np.int(count_time * np.float(divt)), len(longitude))
    )
    Fa_N_down_swap = np.zeros(
        (len(latitude), np.int(count_time * np.float(divt)), len(longitude))
    )
    Fa_N_top_kgpmps_swap = np.swapaxes(Fa_N_top_kgpmps, 0, 1)
    Fa_N_down_kgpmps_swap = np.swapaxes(Fa_N_down_kgpmps, 0, 1)
    for c in range(len(latitude)):
        Fa_N_top_swap[c] = (
            Fa_N_top_kgpmps_swap[c]
            * timestep
            / np.float(divt)
            * 0.5
            * (L_N_gridcell[c] + L_S_gridcell[c])
            / density_water
        )  # [s*m*kg*m^-1*s^-1*kg^-1*m^3]=[m3]
        Fa_N_down_swap[c] = (
            Fa_N_down_kgpmps_swap[c]
            * timestep
            / np.float(divt)
            * 0.5
            * (L_N_gridcell[c] + L_S_gridcell[c])
            / density_water
        )  # [s*m*kg*m^-1*s^-1*kg^-1*m^3]=[m3]

    Fa_N_top_m3 = np.swapaxes(Fa_N_top_swap, 0, 1)
    Fa_N_down_m3 = np.swapaxes(Fa_N_down_swap, 0, 1)

    return Fa_E_top_m3, Fa_E_down_m3, Fa_N_top_m3, Fa_N_down_m3


def get_stablefluxes(
    Fa_E_top,
    Fa_E_down,
    Fa_N_top,
    Fa_N_down,
    timestep,
    divt,
    L_EW_gridcell,
    density_water,
    L_N_gridcell,
    L_S_gridcell,
    latitude,
):

    # find out where the negative fluxes are
    Fa_E_top_posneg = np.ones(np.shape(Fa_E_top))
    Fa_E_top_posneg[Fa_E_top < 0] = -1
    Fa_N_top_posneg = np.ones(np.shape(Fa_E_top))
    Fa_N_top_posneg[Fa_N_top < 0] = -1
    Fa_E_down_posneg = np.ones(np.shape(Fa_E_top))
    Fa_E_down_posneg[Fa_E_down < 0] = -1
    Fa_N_down_posneg = np.ones(np.shape(Fa_E_top))
    Fa_N_down_posneg[Fa_N_down < 0] = -1

    # make everything absolute
    Fa_E_top_abs = np.abs(Fa_E_top)
    Fa_E_down_abs = np.abs(Fa_E_down)
    Fa_N_top_abs = np.abs(Fa_N_top)
    Fa_N_down_abs = np.abs(Fa_N_down)

    # stabilize the outfluxes / influxes
    stab = (
        1.0 / 2.0
    )  # during the reduced timestep the water cannot move further than 1/x * the gridcell,
    # in other words at least x * the reduced timestep is needed to cross a gridcell
    Fa_E_top_stable = np.reshape(
        np.minimum(
            np.reshape(Fa_E_top_abs, (np.size(Fa_E_top_abs))),
            (
                np.reshape(Fa_E_top_abs, (np.size(Fa_E_top_abs)))
                / (
                    np.reshape(Fa_E_top_abs, (np.size(Fa_E_top_abs)))
                    + np.reshape(Fa_N_top_abs, (np.size(Fa_N_top_abs)))
                )
            )
            * stab
            * np.reshape(W_top[:-1, :, :], (np.size(W_top[:-1, :, :]))),
        ),
        (np.int(count_time * np.float(divt)), len(latitude), len(longitude)),
    )
    Fa_N_top_stable = np.reshape(
        np.minimum(
            np.reshape(Fa_N_top_abs, (np.size(Fa_N_top_abs))),
            (
                np.reshape(Fa_N_top_abs, (np.size(Fa_N_top_abs)))
                / (
                    np.reshape(Fa_E_top_abs, (np.size(Fa_E_top_abs)))
                    + np.reshape(Fa_N_top_abs, (np.size(Fa_N_top_abs)))
                )
            )
            * stab
            * np.reshape(W_top[:-1, :, :], (np.size(W_top[:-1, :, :]))),
        ),
        (np.int(count_time * np.float(divt)), len(latitude), len(longitude)),
    )
    Fa_E_down_stable = np.reshape(
        np.minimum(
            np.reshape(Fa_E_down_abs, (np.size(Fa_E_down_abs))),
            (
                np.reshape(Fa_E_down_abs, (np.size(Fa_E_down_abs)))
                / (
                    np.reshape(Fa_E_down_abs, (np.size(Fa_E_down_abs)))
                    + np.reshape(Fa_N_down_abs, (np.size(Fa_N_down_abs)))
                )
            )
            * stab
            * np.reshape(W_down[:-1, :, :], (np.size(W_down[:-1, :, :]))),
        ),
        (np.int(count_time * np.float(divt)), len(latitude), len(longitude)),
    )
    Fa_N_down_stable = np.reshape(
        np.minimum(
            np.reshape(Fa_N_down_abs, (np.size(Fa_N_down_abs))),
            (
                np.reshape(Fa_N_down_abs, (np.size(Fa_N_down_abs)))
                / (
                    np.reshape(Fa_E_down_abs, (np.size(Fa_E_down_abs)))
                    + np.reshape(Fa_N_down_abs, (np.size(Fa_N_down_abs)))
                )
            )
            * stab
            * np.reshape(W_down[:-1, :, :], (np.size(W_down[:-1, :, :]))),
        ),
        (np.int(count_time * np.float(divt)), len(latitude), len(longitude)),
    )

    # get rid of the nan values
    Fa_E_top_stable[np.isnan(Fa_E_top_stable)] = 0
    Fa_N_top_stable[np.isnan(Fa_N_top_stable)] = 0
    Fa_E_down_stable[np.isnan(Fa_E_down_stable)] = 0
    Fa_N_down_stable[np.isnan(Fa_N_down_stable)] = 0

    # redefine
    Fa_E_top = Fa_E_top_stable * Fa_E_top_posneg
    Fa_N_top = Fa_N_top_stable * Fa_N_top_posneg
    Fa_E_down = Fa_E_down_stable * Fa_E_down_posneg
    Fa_N_down = Fa_N_down_stable * Fa_N_down_posneg

    return Fa_E_top, Fa_E_down, Fa_N_top, Fa_N_down


# Code
def getFa_Vert(
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
):

    # total moisture in the column
    W = W_top + W_down

    # define the horizontal fluxes over the boundaries
    # fluxes over the eastern boundary
    Fa_E_top_boundary = np.zeros(np.shape(Fa_E_top))
    Fa_E_top_boundary[:, :, :-1] = 0.5 * (Fa_E_top[:, :, :-1] + Fa_E_top[:, :, 1:])
    Fa_E_down_boundary = np.zeros(np.shape(Fa_E_down))
    Fa_E_down_boundary[:, :, :-1] = 0.5 * (Fa_E_down[:, :, :-1] + Fa_E_down[:, :, 1:])

    # find out where the positive and negative fluxes are
    Fa_E_top_pos = np.ones(np.shape(Fa_E_top))
    Fa_E_down_pos = np.ones(np.shape(Fa_E_down))
    Fa_E_top_pos[Fa_E_top_boundary < 0] = 0
    Fa_E_down_pos[Fa_E_down_boundary < 0] = 0
    Fa_E_top_neg = Fa_E_top_pos - 1
    Fa_E_down_neg = Fa_E_down_pos - 1

    # separate directions west-east (all positive numbers)
    Fa_E_top_WE = Fa_E_top_boundary * Fa_E_top_pos
    Fa_E_top_EW = Fa_E_top_boundary * Fa_E_top_neg
    Fa_E_down_WE = Fa_E_down_boundary * Fa_E_down_pos
    Fa_E_down_EW = Fa_E_down_boundary * Fa_E_down_neg

    # fluxes over the western boundary
    Fa_W_top_WE = np.nan * np.zeros(np.shape(P))
    Fa_W_top_WE[:, :, 1:] = Fa_E_top_WE[:, :, :-1]
    Fa_W_top_WE[:, :, 0] = Fa_E_top_WE[:, :, -1]
    Fa_W_top_EW = np.nan * np.zeros(np.shape(P))
    Fa_W_top_EW[:, :, 1:] = Fa_E_top_EW[:, :, :-1]
    Fa_W_top_EW[:, :, 0] = Fa_E_top_EW[:, :, -1]
    Fa_W_down_WE = np.nan * np.zeros(np.shape(P))
    Fa_W_down_WE[:, :, 1:] = Fa_E_down_WE[:, :, :-1]
    Fa_W_down_WE[:, :, 0] = Fa_E_down_WE[:, :, -1]
    Fa_W_down_EW = np.nan * np.zeros(np.shape(P))
    Fa_W_down_EW[:, :, 1:] = Fa_E_down_EW[:, :, :-1]
    Fa_W_down_EW[:, :, 0] = Fa_E_down_EW[:, :, -1]

    # fluxes over the northern boundary
    Fa_N_top_boundary = np.nan * np.zeros(np.shape(Fa_N_top))
    Fa_N_top_boundary[:, 1:, :] = 0.5 * (Fa_N_top[:, :-1, :] + Fa_N_top[:, 1:, :])
    Fa_N_down_boundary = np.nan * np.zeros(np.shape(Fa_N_down))
    Fa_N_down_boundary[:, 1:, :] = 0.5 * (Fa_N_down[:, :-1, :] + Fa_N_down[:, 1:, :])

    # find out where the positive and negative fluxes are
    Fa_N_top_pos = np.ones(np.shape(Fa_N_top))
    Fa_N_down_pos = np.ones(np.shape(Fa_N_down))
    Fa_N_top_pos[Fa_N_top_boundary < 0] = 0
    Fa_N_down_pos[Fa_N_down_boundary < 0] = 0
    Fa_N_top_neg = Fa_N_top_pos - 1
    Fa_N_down_neg = Fa_N_down_pos - 1

    # separate directions south-north (all positive numbers)
    Fa_N_top_SN = Fa_N_top_boundary * Fa_N_top_pos
    Fa_N_top_NS = Fa_N_top_boundary * Fa_N_top_neg
    Fa_N_down_SN = Fa_N_down_boundary * Fa_N_down_pos
    Fa_N_down_NS = Fa_N_down_boundary * Fa_N_down_neg

    # fluxes over the southern boundary
    Fa_S_top_SN = np.nan * np.zeros(np.shape(P))
    Fa_S_top_SN[:, :-1, :] = Fa_N_top_SN[:, 1:, :]
    Fa_S_top_NS = np.nan * np.zeros(np.shape(P))
    Fa_S_top_NS[:, :-1, :] = Fa_N_top_NS[:, 1:, :]
    Fa_S_down_SN = np.nan * np.zeros(np.shape(P))
    Fa_S_down_SN[:, :-1, :] = Fa_N_down_SN[:, 1:, :]
    Fa_S_down_NS = np.nan * np.zeros(np.shape(P))
    Fa_S_down_NS[:, :-1, :] = Fa_N_down_NS[:, 1:, :]

    # check the water balance
    Sa_after_Fa_down = np.zeros([1, len(latitude), len(longitude)])
    Sa_after_Fa_top = np.zeros([1, len(latitude), len(longitude)])
    Sa_after_all_down = np.zeros([1, len(latitude), len(longitude)])
    Sa_after_all_top = np.zeros([1, len(latitude), len(longitude)])
    residual_down = np.zeros(np.shape(P))  # residual factor [m3]
    residual_top = np.zeros(np.shape(P))  # residual factor [m3]

    for t in range(np.int(count_time * divt)):
        # down: calculate with moisture fluxes:
        Sa_after_Fa_down[0, 1:-1, :] = (
            W_down[t, 1:-1, :]
            - Fa_E_down_WE[t, 1:-1, :]
            + Fa_E_down_EW[t, 1:-1, :]
            + Fa_W_down_WE[t, 1:-1, :]
            - Fa_W_down_EW[t, 1:-1, :]
            - Fa_N_down_SN[t, 1:-1, :]
            + Fa_N_down_NS[t, 1:-1, :]
            + Fa_S_down_SN[t, 1:-1, :]
            - Fa_S_down_NS[t, 1:-1, :]
        )

        # top: calculate with moisture fluxes:
        Sa_after_Fa_top[0, 1:-1, :] = (
            W_top[t, 1:-1, :]
            - Fa_E_top_WE[t, 1:-1, :]
            + Fa_E_top_EW[t, 1:-1, :]
            + Fa_W_top_WE[t, 1:-1, :]
            - Fa_W_top_EW[t, 1:-1, :]
            - Fa_N_top_SN[t, 1:-1, :]
            + Fa_N_top_NS[t, 1:-1, :]
            + Fa_S_top_SN[t, 1:-1, :]
            - Fa_S_top_NS[t, 1:-1, :]
        )

        # down: substract precipitation and add evaporation
        Sa_after_all_down[0, 1:-1, :] = (
            Sa_after_Fa_down[0, 1:-1, :]
            - P[t, 1:-1, :] * (W_down[t, 1:-1, :] / W[t, 1:-1, :])
            + E[t, 1:-1, :]
        )

        # top: substract precipitation
        Sa_after_all_top[0, 1:-1, :] = Sa_after_Fa_top[0, 1:-1, :] - P[t, 1:-1, :] * (
            W_top[t, 1:-1, :] / W[t, 1:-1, :]
        )

        # down: calculate the residual
        residual_down[t, 1:-1, :] = (
            W_down[t + 1, 1:-1, :] - Sa_after_all_down[0, 1:-1, :]
        )

        # top: calculate the residual
        residual_top[t, 1:-1, :] = W_top[t + 1, 1:-1, :] - Sa_after_all_top[0, 1:-1, :]

    # compute the resulting vertical moisture flux
    Fa_Vert_raw = (
        W_down[1:, :, :] / W[1:, :, :] * (residual_down + residual_top) - residual_down
    )  # the vertical velocity so that the new residual_down/W_down =  residual_top/W_top (positive downward)

    # find out where the negative vertical flux is
    Fa_Vert_posneg = np.ones(np.shape(Fa_Vert_raw))
    Fa_Vert_posneg[Fa_Vert_raw < 0] = -1

    # make the vertical flux absolute
    Fa_Vert_abs = np.abs(Fa_Vert_raw)

    # stabilize the outfluxes / influxes
    stab = (
        1.0 / 4.0
    )  # during the reduced timestep the vertical flux can maximally empty/fill 1/x of the top or down storage

    Fa_Vert_stable = np.reshape(
        np.minimum(
            np.reshape(Fa_Vert_abs, (np.size(Fa_Vert_abs))),
            np.minimum(
                stab * np.reshape(W_top[1:, :, :], (np.size(W_top[1:, :, :]))),
                stab * np.reshape(W_down[1:, :, :], (np.size(W_down[1:, :, :]))),
            ),
        ),
        (np.int(count_time * np.float(divt)), len(latitude), len(longitude)),
    )

    # redefine the vertical flux
    Fa_Vert = Fa_Vert_stable * Fa_Vert_posneg

    return Fa_Vert_raw, Fa_Vert, residual_down, residual_top


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
    #start = timer()
    print(("0 = " + str(timer())))
    #    #1 integrate specific humidity to get the (total) column water (vapor) and calculate horizontal moisture fluxes
    cwv, W_top, W_down, Fa_E_top, Fa_N_top, Fa_E_down, Fa_N_down = getWandFluxes(
        latnrs,
        lonnrs,
        date,
        density_water,
        latitude,
        longitude,
        g,
        A_gridcell,
    )
    print(("1,2,3 = " + str(timer())))

    # 4 evaporation and precipitation
    # daily loop
    E, P = getEP(
        latnrs,
        lonnrs,
        yearnumber,
        monthnumber,
        begin_time,
        count_time,
        latitude,
        longitude,
        A_gridcell,
    ) # uses hourly data
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
