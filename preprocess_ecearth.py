import calendar
import datetime as dt
import os

# Import libraries
import numpy as np
import pandas as pd
import xarray as xr
import scipy.io as sio
import yaml

from scipy.interpolate import interp1d

from getconstants_pressure_ECEarth import getconstants_pressure_ECEarth

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
datelist = pd.date_range(
    start=config['start_date'],
    end=config['end_date'],
    freq="d",
    inclusive="left"
    )


def get_input_data(variable, year, month):
    """Get input data for variable."""
    filename = f"{name_of_run}{variable}_{year}{month:02d}_NH.nc"
    filepath = os.path.join(config["input_folder"], filename)
    return xr.open_dataset(filepath)


def get_input_data_next_month(variable, year, month):
    """Get input data for next month."""
    if month == 12:
        return get_input_data(variable, year + 1, 1)
    else:
        return get_input_data(variable, year, month + 1)


def get_output_path(year, month, a):
    """Get data for output file."""
    filename = f"{year}-{month:02d}-{a:02d}fluxes_storages.mat"
    save_path = os.path.join(config["interdata_folder"], filename)
    return save_path


# Determine the fluxes and states
# In this defintion the vertical spline interpolation is performed to determine the moisture fluxes for the two layers at each grid cell
def getWandFluxes(
    latnrs,
    lonnrs,
    final_time,
    a,
    year,
    month,
    begin_time,
    count_time,
    density_water,
    latitude,
    longitude,
    g,
    A_gridcell,
):

    times_6hourly = slice(0, 5)
    times_3hourly = slice(0, 10, 2)

    q = get_input_data("Q", year, month).Q.isel(time=times_6hourly, lat=latnrs, lon=lonnrs)
    u = get_input_data("U", year, month).U.isel(time=times_6hourly, lat=latnrs, lon=lonnrs)
    v = get_input_data("V", year, month).V.isel(time=times_6hourly, lat=latnrs, lon=lonnrs)
    q2m = get_input_data("Q2M", year, month).Q2M.isel(time=times_3hourly, lat=latnrs, lon=lonnrs)
    lnsp = get_input_data("LNSP", year, month).LNSP.isel(time=times_3hourly, lat=latnrs, lon=lonnrs)
    u10 = get_input_data("U10", year, month).U10M.isel(time=times_6hourly, lat=latnrs, lon=lonnrs)
    v10 = get_input_data("V10", year, month).V10M.isel(time=times_6hourly, lat=latnrs, lon=lonnrs)

    if a == final_time:
        q_next = get_input_data_next_month("Q", year, month).Q.isel(time=0, lat=latnrs, lon=lonnrs)
        u_next = get_input_data_next_month("U", year, month).U.isel(time=0, lat=latnrs, lon=lonnrs)
        v_next = get_input_data_next_month("V", year, month).V.isel(time=0, lat=latnrs, lon=lonnrs)
        q2m_next = get_input_data_next_month("Q2M", year, month).Q2M.isel(time=0, lat=latnrs, lon=lonnrs)
        lnsp_next = get_input_data_next_month("LNSP", year, month).LNSP.isel(time=0, lat=latnrs, lon=lonnrs)
        u10_next = get_input_data_next_month("U10", year, month).U10M.isel(time=0, lat=latnrs, lon=lonnrs)
        v10_next = get_input_data_next_month("V10", year, month).V10M.isel(time=0, lat=latnrs, lon=lonnrs)

        q = xr.concat([q, q_next], dim='time')
        u = xr.concat([u, u_next], dim='time')
        v = xr.concat([v, v_next], dim='time')
        q2m = xr.concat([q2m, q2m_next], dim='time')
        lnsp = xr.concat([lnsp, lnsp_next], dim='time')
        u10 = xr.concat([u10, u10_next], dim='time')
        v10 = xr.concat([v10, v10_next], dim='time')

    sp = np.exp(lnsp)  # log(sp) --> sp (Pa)

    # For now extract bare numpy arrays from the data
    levelist = u.lev.values
    q = q.values
    u = u.values
    v = v.values
    q2m = q2m.values
    sp = sp.values.squeeze() # drop len(1) dimension 'lev'
    u10 = u10.values
    v10 = v10.values

    n_levels = 40  # from five pressure levels the vertical data is interpolated to 40 levels

    ##Imme: location of boundary is hereby hard defined at model level 47 which corresponds with about
    P_boundary = 0.72878581 * sp + 7438.803223
    dp = (sp - 20000.0) / (n_levels - 1)
    time = q.shape[0]

    p_maxmin = np.zeros((time, n_levels + 2, len(latitude), len(longitude)))
    p_maxmin[:, 1:-1, :, :] = (sp[:, np.newaxis, :, :] - dp[:, np.newaxis, :, :] * np.arange(0, n_levels)[np.newaxis, :, np.newaxis, np.newaxis])

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
        "after p_maxmin and now",
        dt.datetime.now().time(),
    )

    # Add surface and atmosphere together for u,q,v
    p = np.zeros((time, levelist.size + 2, len(latitude), len(longitude)))
    p[:, 1:-1, :, :] = levelist[np.newaxis, :, np.newaxis, np.newaxis]
    p[:, 0, :, :] = sp

    u_total = np.zeros((time, levelist.size + 2, len(latitude), len(longitude)))
    u_total[:, 1:-1, :, :] = u
    u_total[:, 0, :, :] = u10
    u_total[:, -1, :, :] = u[:, -1, :, :]

    v_total = np.zeros((time, levelist.size + 2, len(latitude), len(longitude)))
    v_total[:, 1:-1, :, :] = v
    v_total[:, 0, :, :] = v10
    v_total[:, -1, :, :] = v[:, -1, :, :]

    q_total = np.zeros((time, levelist.size + 2, len(latitude), len(longitude)))
    q_total[:, 1:-1, :, :] = q
    q_total[:, 0, :, :] = q2m

    mask = np.ones(u_total.shape, dtype=np.bool)
    mask[:, 1:-1, :, :] = levelist[np.newaxis, :, np.newaxis, np.newaxis] < (
        sp[:, np.newaxis, :, :] - 1000.0
    )  # Pa

    u_masked = np.ma.masked_array(u_total, mask=~mask)
    v_masked = np.ma.masked_array(v_total, mask=~mask)
    q_masked = np.ma.masked_array(q_total, mask=~mask)
    p_masked = np.ma.masked_array(p, mask=~mask)

    del (u_total, v_total, q_total, p, u, v, q, u10, v10, q2m, sp)

    print("Data is loaded", dt.datetime.now().time())
    print("before interpolation loop", dt.datetime.now().time())

    uq_maxmin = np.zeros((time, n_levels + 2, len(latitude), len(longitude)))
    vq_maxmin = np.zeros((time, n_levels + 2, len(latitude), len(longitude)))
    q_maxmin = np.zeros((time, n_levels + 2, len(latitude), len(longitude)))

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


def getEP(
    latnrs, lonnrs, year, month, begin_time, count_time, latitude, longitude, A_gridcell
):
    # 3-hourly data so 8 steps per day

    # (accumulated after the forecast at 00.00 and 12.00 by steps of 3 hours in time
    evaporation = get_input_data("EVAP", year, month).variables["E"][
        begin_time * 2 : (begin_time * 2 + count_time * 2), latnrs, lonnrs
    ]  # m
    precipitation = get_input_data("TP", year, month).variables["TP"][
        begin_time * 2 : (begin_time * 2 + count_time * 2), latnrs, lonnrs
    ]  # m

    # delete and transfer negative values, change sign convention to all positive
    precipitation = np.reshape(
        np.maximum(
            np.reshape(precipitation, (np.size(precipitation)))
            + np.maximum(np.reshape(evaporation, (np.size(evaporation))), 0.0),
            0.0,
        ),
        (np.int(count_time * 2), len(latitude), len(longitude)),
    )
    evaporation = np.reshape(
        np.abs(np.minimum(np.reshape(evaporation, (np.size(evaporation))), 0.0)),
        (np.int(count_time * 2), len(latitude), len(longitude)),
    )

    # calculate volumes
    A_gridcell2D = np.tile(A_gridcell, [1, len(longitude)])
    A_gridcell_1_2D = np.reshape(A_gridcell2D, [1, len(latitude), len(longitude)])
    A_gridcell_max3D = np.tile(A_gridcell_1_2D, [count_time * 2, 1, 1])

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

start1 = dt.datetime.now()
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
) = getconstants_pressure_ECEarth(latnrs, lonnrs, config["land_sea_mask"])

def time_since(start_time):
    """Return elapsed time in a formatted string."""
    elapsed_time = dt.datetime.now() - start_time
    return f"{minutes:.0f}:{seconds:.02f}"

for date in datelist:
    start = dt.datetime.now()

    day = date.day
    year = date.year
    month = date.month

    begin_time = (day - 1) * count_time  # Python starts at 0
    final_time = calendar.monthrange(year, month)[1]

    print(f"Preprocessing data for {date}, {begin_time=}, {final_time=}")

    # 1 integrate specific humidity to get the (total) column water (vapor) and calculate horizontal moisture fluxes
    cwv, W_top, W_down, Fa_E_top, Fa_N_top, Fa_E_down, Fa_N_down = getWandFluxes(
        latnrs,
        lonnrs,
        final_time,
        day,
        year,
        month,
        begin_time,
        count_time,
        density_water,
        latitude,
        longitude,
        g,
        A_gridcell,
    )
    print(f"Step 1, 2, 3 finished, elapsed time since start: {dt.datetime.now() - start}")

    # 4 evaporation and precipitation
    E, P = getEP(
        latnrs,
        lonnrs,
        year,
        month,
        begin_time,
        count_time,
        latitude,
        longitude,
        A_gridcell,
    )
    print(f"Step 4 finished, elapsed time since start: {dt.datetime.now() - start}")

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
    print(f"Step 5 finished, elapsed time since start: {dt.datetime.now() - start}")

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
    print(f"Step 6a finished, elapsed time since start: {dt.datetime.now() - start}")

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
    print(f"Step 6b finished, elapsed time since start: {dt.datetime.now() - start}")

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
    print(f"Step 7 finished, elapsed time since start: {dt.datetime.now() - start}")

    sio.savemat(
        get_output_path(year, month, day),
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

    end = dt.datetime.now()
    print(f"Runtime fluxes_and_storages for {day=}, {year=} is {end-start}")

end1 = dt.datetime.now()
print("The total runtime is {end1 - start1}")
