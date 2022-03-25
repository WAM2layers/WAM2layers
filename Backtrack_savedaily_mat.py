# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 13:24:45 2016

@author: Ent00002
"""

"""
Created on Mon Feb 18 15:30:43 2019

@author: bened003
"""

# This Backtrack script is based on the Con_E_Recyc_Masterscript.py from the WAM-2layers model from Ruud van der Ent but adapted to fit the purposes of the study by Imme Benedict

# We have removed the water_lost statement to conserve water mass
# We have implemented a datelist function so the model can run for multiple years without having problems with leap years
# Timetracking was not implemented in this code

# We save daily data instead of data at every timestep to reduce storage

import calendar
import datetime as dt
import os
from datetime import timedelta
from timeit import default_timer as timer

# Import libraries
import numpy as np
import scipy.io as sio

from getconstants_pressure_ECEarth import getconstants_pressure_ECEarth


# to create datelist
def get_times_daily(startdate, enddate):
    """generate a dictionary with date/times"""
    numdays = enddate - startdate
    dateList = []
    for x in range(0, numdays.days + 1):
        dateList.append(startdate + dt.timedelta(days=x))
    return dateList


# BEGIN OF INPUT1 (FILL THIS IN)

months = np.arange(1, 13)  # for full year, enter np.arange(1,13)
months_length_leap = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
months_length_nonleap = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
years = np.arange(
    2002, 2007
)  # fill in the years # If I fill in more than one year than I need to set the months to 12 (so according to python counting to 13)

# create datelist
if int(calendar.isleap(years[-1])) == 0:  # no leap year
    datelist = get_times_daily(
        dt.date(years[0], months[0], 1),
        dt.date(years[-1], months[-1], months_length_nonleap[months[-1] - 1]),
    )[::-1]
else:  # leap
    datelist = get_times_daily(
        dt.date(years[0], months[0], 1),
        dt.date(years[-1], months[-1], months_length_leap[months[-1] - 1]),
    )[
        ::-1
    ]  # [::-1] to run the datelist the other way around

divt = 60  # division of the timestep, 96 means a calculation timestep of 24/96 = 0.25 hours (numerical stability purposes)
count_time = 4  # number of indices to get data from (for daily data this means everytime one day)

# Manage the extent of your dataset (FILL THIS IN)
latnrs = np.arange(0, 267)  # minimal domain
lonnrs = np.arange(0, 444)

isglobal = 0  # fill in 1 for global computations (i.e. Earth round), fill in 0 for a local domain with boundaries

# obtain the constants
lsm_data_ECEarth_T799 = "landseamask_ECearth_T799.nc"  # insert landseamask here

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
) = getconstants_pressure_ECEarth(latnrs, lonnrs, lsm_data_ECEarth_T799)

# BEGIN OF INPUT 2 (FILL THIS IN)
Region = np.load("mask_Miss_ECEarth.npy")  # region to perform the tracking for
Kvf = 3  # vertical dispersion factor (advection only is 0, dispersion the same size of the advective flux is 1, for stability don't make this more than 3)
timetracking = 0  # 0 for not tracking time and 1 for tracking time
veryfirstrun = 1  # type '1' if no run has been done before from which can be continued, otherwise type '0'

# END OF INPUT
# Datapaths (FILL THIS IN)

interdata_folder = "Interdata_ECEarth/PresentMember2_correct/"  # must be an existing folder # insert Interdata folder here

# Check if interdata folder exists:
assert os.path.isdir(
    interdata_folder
), "Please create the interdata_folder before running the script"
# Check if sub interdata folder exists otherwise create it:
sub_interdata_folder = os.path.join(interdata_folder, "Regional_backward_daily")

if os.path.isdir(sub_interdata_folder):
    pass
else:
    os.makedirs(sub_interdata_folder)


def data_path_ea(year, month, day):
    save_empty_arrays_track = os.path.join(
        sub_interdata_folder,
        str(year) + "-" + str(month).zfill(2) + "-" + str(day).zfill(2) + "Sa_track",
    )
    save_empty_arrays_time = os.path.join(
        sub_interdata_folder,
        str(year) + "-" + str(month).zfill(2) + "-" + str(day).zfill(2) + "Sa_time",
    )
    return save_empty_arrays_track, save_empty_arrays_time


def data_path(previous_data_to_load, yearnumber, month, a):
    load_Sa_track = os.path.join(
        sub_interdata_folder,
        str(previous_data_to_load.year)
        + "-"
        + str(previous_data_to_load.month).zfill(2)
        + "-"
        + str(previous_data_to_load.day).zfill(2)
        + "Sa_track.npz",
    )

    load_fluxes_and_storages = os.path.join(
        interdata_folder,
        str(yearnumber)
        + "-"
        + str(month).zfill(2)
        + "-"
        + str(a).zfill(2)
        + "fluxes_storages.mat",
    )

    load_Sa_time = os.path.join(
        sub_interdata_folder,
        str(previous_data_to_load.year)
        + "-"
        + str(previous_data_to_load.month).zfill(2)
        + "-"
        + str(previous_data_to_load.day).zfill(2)
        + "Sa_time.npz",
    )

    save_path_track = os.path.join(
        sub_interdata_folder,
        str(yearnumber)
        + "-"
        + str(month).zfill(2)
        + "-"
        + str(a).zfill(2)
        + "Sa_track",
    )
    save_path_time = os.path.join(
        sub_interdata_folder,
        str(yearnumber) + "-" + str(month).zfill(2) + "-" + str(a).zfill(2) + "Sa_time",
    )
    return (
        load_Sa_track,
        load_fluxes_and_storages,
        load_Sa_time,
        save_path_track,
        save_path_time,
    )


# Code


def get_Sa_track_backward(
    latitude,
    longitude,
    count_time,
    divt,
    Kvf,
    Region,
    Fa_E_top,
    Fa_N_top,
    Fa_E_down,
    Fa_N_down,
    Fa_Vert,
    E,
    P,
    W_top,
    W_down,
    Sa_track_top_last,
    Sa_track_down_last,
):

    # make P_region matrix
    Region3D = np.tile(
        np.reshape(Region, [1, len(latitude), len(longitude)]), [len(P[:, 0, 0]), 1, 1]
    )
    P_region = Region3D * P

    # Total moisture in the column
    W = W_top + W_down

    # separate the direction of the vertical flux and make it absolute
    Fa_upward = np.zeros(np.shape(Fa_Vert))
    Fa_upward[Fa_Vert <= 0] = Fa_Vert[
        Fa_Vert <= 0
    ]  # in 4th timestep: __main__:1: RuntimeWarning: invalid value encountered in less_equal # not a problem
    Fa_downward = np.zeros(np.shape(Fa_Vert))
    Fa_downward[Fa_Vert >= 0] = Fa_Vert[Fa_Vert >= 0]
    Fa_upward = np.abs(Fa_upward)

    # include the vertical dispersion
    if Kvf == 0:
        pass
        # do nothing
    else:
        Fa_upward = (1.0 + Kvf) * Fa_upward
        Fa_upward[Fa_Vert >= 0] = Fa_Vert[Fa_Vert >= 0] * Kvf
        Fa_downward = (1.0 + Kvf) * Fa_downward
        Fa_downward[Fa_Vert <= 0] = np.abs(Fa_Vert[Fa_Vert <= 0]) * Kvf

    # define the horizontal fluxes over the boundaries
    # fluxes over the eastern boundary
    Fa_E_top_boundary = np.zeros(np.shape(Fa_E_top))
    Fa_E_top_boundary[:, :, :-1] = 0.5 * (Fa_E_top[:, :, :-1] + Fa_E_top[:, :, 1:])
    if isglobal == 1:
        Fa_E_top_boundary[:, :, -1] = 0.5 * (Fa_E_top[:, :, -1] + Fa_E_top[:, :, 0])
    Fa_E_down_boundary = np.zeros(np.shape(Fa_E_down))
    Fa_E_down_boundary[:, :, :-1] = 0.5 * (Fa_E_down[:, :, :-1] + Fa_E_down[:, :, 1:])
    if isglobal == 1:
        Fa_E_down_boundary[:, :, -1] = 0.5 * (Fa_E_down[:, :, -1] + Fa_E_down[:, :, 0])

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
    # Fa_N_top_boundary = np.nan*np.zeros(np.shape(Fa_N_top)); #Imme: why do you multiply here your zeros with nans ?!?!?!? you do not do that for Fa_E_top_boundary
    # Adapted by Imme
    Fa_N_top_boundary = np.zeros(np.shape(Fa_N_top))
    # Imme: why do you multiply here your zeros with nans ?!?!?!? you do not do that for Fa_E_top_boundary
    Fa_N_top_boundary[:, 1:, :] = 0.5 * (Fa_N_top[:, :-1, :] + Fa_N_top[:, 1:, :])
    # Fa_N_down_boundary = np.nan*np.zeros(np.shape(Fa_N_down));
    # Adapted by Imme
    Fa_N_down_boundary = np.zeros(np.shape(Fa_N_down))
    # als je niet met np.nan multiplied krijg je in de volgende alinea geen invalid value encountered in less
    # verandert er verder nog wat!??!
    Fa_N_down_boundary[:, 1:, :] = 0.5 * (Fa_N_down[:, :-1, :] + Fa_N_down[:, 1:, :])

    # find out where the positive and negative fluxes are
    Fa_N_top_pos = np.ones(np.shape(Fa_N_top))
    Fa_N_down_pos = np.ones(np.shape(Fa_N_down))
    Fa_N_top_pos[
        Fa_N_top_boundary < 0
    ] = 0  # Invalid value encountered in less, omdat er nan in Fa_N_top_boundary staan
    Fa_N_down_pos[
        Fa_N_down_boundary < 0
    ] = 0  # Invalid value encountered in less, omdat er nan in Fa_N_top_boundary staan
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

    # defining size of output
    Sa_track_down = np.zeros(np.shape(W_down))
    Sa_track_top = np.zeros(np.shape(W_top))

    # assign begin values of output == last (but first index) values of the previous time slot
    Sa_track_down[-1, :, :] = Sa_track_down_last
    Sa_track_top[-1, :, :] = Sa_track_top_last

    # defining sizes of tracked moisture
    Sa_track_after_Fa_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_track_after_Fa_P_E_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_track_E_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_track_W_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_track_N_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_track_S_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_track_after_Fa_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_track_after_Fa_P_E_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_track_E_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_track_W_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_track_N_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_track_S_top = np.zeros(np.shape(Sa_track_top_last))

    # define sizes of total moisture
    Sa_E_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_W_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_N_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_S_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_E_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_W_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_N_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_S_top = np.zeros(np.shape(Sa_track_top_last))

    # define variables that find out what happens to the water
    north_loss = np.zeros((np.int(count_time * divt), 1, len(longitude)))
    south_loss = np.zeros((np.int(count_time * divt), 1, len(longitude)))
    east_loss = np.zeros((np.int(count_time * divt), 1, len(latitude)))
    west_loss = np.zeros((np.int(count_time * divt), 1, len(latitude)))
    down_to_top = np.zeros(np.shape(P))
    top_to_down = np.zeros(np.shape(P))
    water_lost = np.zeros(np.shape(P))
    water_lost_down = np.zeros(np.shape(P))
    water_lost_top = np.zeros(np.shape(P))

    # Sa calculation backward in time
    for t in np.arange(np.int(count_time * divt), 0, -1):
        # down: define values of total moisture
        Sa_E_down[0, :, :-1] = W_down[
            t, :, 1:
        ]  # Atmospheric storage of the cell to the east [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors
        # Sa_E_down[0,:,-1] = W_down[t,:,0] # Atmospheric storage of the cell to the east [m3]
        Sa_W_down[0, :, 1:] = W_down[
            t, :, :-1
        ]  # Atmospheric storage of the cell to the west [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors
        # Sa_W_down[0,:,0] = W_down[t,:,-1] # Atmospheric storage of the cell to the west [m3]
        Sa_N_down[0, 1:, :] = W_down[
            t, 0:-1, :
        ]  # Atmospheric storage of the cell to the north [m3]
        Sa_S_down[0, :-1, :] = W_down[
            t, 1:, :
        ]  # Atmospheric storage of the cell to the south [m3]

        # top: define values of total moisture
        Sa_E_top[0, :, :-1] = W_top[
            t, :, 1:
        ]  # Atmospheric storage of the cell to the east [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors
        # Sa_E_top[0,:,-1] = W_top[t,:,0] # Atmospheric storage of the cell to the east [m3]
        Sa_W_top[0, :, 1:] = W_top[
            t, :, :-1
        ]  # Atmospheric storage of the cell to the west [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors
        # Sa_W_top[0,:,0] = W_top[t,:,-1] # Atmospheric storage of the cell to the west [m3]
        Sa_N_top[0, 1:, :] = W_top[
            t, :-1, :
        ]  # Atmospheric storage of the cell to the north [m3]
        Sa_S_top[0, :-1, :] = W_top[
            t, 1:, :
        ]  # Atmospheric storage of the cell to the south [m3]

        # down: define values of tracked moisture of neighbouring grid cells
        Sa_track_E_down[0, :, :-1] = Sa_track_down[
            t, :, 1:
        ]  # Atmospheric tracked storage of the cell to the east [m3]
        if isglobal == 1:
            Sa_track_E_down[0, :, -1] = Sa_track_down[
                t, :, 0
            ]  # Atmospheric tracked storage of the cell to the east [m3]
        Sa_track_W_down[0, :, 1:] = Sa_track_down[
            t, :, :-1
        ]  # Atmospheric storage of the cell to the west [m3]
        if isglobal == 1:
            Sa_track_W_down[0, :, 0] = Sa_track_down[
                t, :, -1
            ]  # Atmospheric storage of the cell to the west [m3]
        Sa_track_N_down[0, 1:, :] = Sa_track_down[
            t, :-1, :
        ]  # Atmospheric storage of the cell to the north [m3]
        Sa_track_S_down[0, :-1, :] = Sa_track_down[
            t, 1:, :
        ]  # Atmospheric storage of the cell to the south [m3]

        # down: calculate with moisture fluxes
        Sa_track_after_Fa_down[0, 1:-1, 1:-1] = (
            Sa_track_down[t, 1:-1, 1:-1]
            + Fa_E_down_WE[t - 1, 1:-1, 1:-1]
            * (Sa_track_E_down[0, 1:-1, 1:-1] / Sa_E_down[0, 1:-1, 1:-1])
            - Fa_E_down_EW[t - 1, 1:-1, 1:-1]
            * (Sa_track_down[t, 1:-1, 1:-1] / W_down[t, 1:-1, 1:-1])
            - Fa_W_down_WE[t - 1, 1:-1, 1:-1]
            * (Sa_track_down[t, 1:-1, 1:-1] / W_down[t, 1:-1, 1:-1])
            + Fa_W_down_EW[t - 1, 1:-1, 1:-1]
            * (Sa_track_W_down[0, 1:-1, 1:-1] / Sa_W_down[0, 1:-1, 1:-1])
            + Fa_N_down_SN[t - 1, 1:-1, 1:-1]
            * (Sa_track_N_down[0, 1:-1, 1:-1] / Sa_N_down[0, 1:-1, 1:-1])
            - Fa_N_down_NS[t - 1, 1:-1, 1:-1]
            * (Sa_track_down[t, 1:-1, 1:-1] / W_down[t, 1:-1, 1:-1])
            - Fa_S_down_SN[t - 1, 1:-1, 1:-1]
            * (Sa_track_down[t, 1:-1, 1:-1] / W_down[t, 1:-1, 1:-1])
            + Fa_S_down_NS[t - 1, 1:-1, 1:-1]
            * (Sa_track_S_down[0, 1:-1, 1:-1] / Sa_S_down[0, 1:-1, 1:-1])
            - Fa_downward[t - 1, 1:-1, 1:-1]
            * (Sa_track_down[t, 1:-1, 1:-1] / W_down[t, 1:-1, 1:-1])
            + Fa_upward[t - 1, 1:-1, 1:-1]
            * (Sa_track_top[t, 1:-1, 1:-1] / W_top[t, 1:-1, 1:-1])
        )
        # sometimes you get a Runtimewarning: invalid value encountered in less
        # RuntimeWarning: invalid value encountered in divide
        # ik kan me voorstellen dat die error er komt als er een nan in Sa_track_down zit en dat je daar dan door moet gaan delen..

        # top: define values of tracked moisture of neighbouring grid cells
        Sa_track_E_top[0, :, :-1] = Sa_track_top[
            t, :, 1:
        ]  # Atmospheric tracked storage of the cell to the east [m3]
        if isglobal == 1:
            Sa_track_E_top[0, :, -1] = Sa_track_top[
                t, :, 0
            ]  # Atmospheric tracked storage of the cell to the east [m3]
        Sa_track_W_top[0, :, 1:] = Sa_track_top[
            t, :, :-1
        ]  # Atmospheric tracked storage of the cell to the west [m3]
        if isglobal == 1:
            Sa_track_W_top[0, :, 0] = Sa_track_top[
                t, :, -1
            ]  # Atmospheric tracked storage of the cell to the west [m3]
        Sa_track_N_top[0, 1:, :] = Sa_track_top[
            t, :-1, :
        ]  # Atmospheric tracked storage of the cell to the north [m3]
        Sa_track_S_top[0, :-1, :] = Sa_track_top[
            t, 1:, :
        ]  # Atmospheric tracked storage of the cell to the south [m3]

        # top: calculate with moisture fluxes
        Sa_track_after_Fa_top[0, 1:-1, 1:-1] = (
            Sa_track_top[t, 1:-1, 1:-1]
            + Fa_E_top_WE[t - 1, 1:-1, 1:-1]
            * (Sa_track_E_top[0, 1:-1, 1:-1] / Sa_E_top[0, 1:-1, 1:-1])
            - Fa_E_top_EW[t - 1, 1:-1, 1:-1]
            * (Sa_track_top[t, 1:-1, 1:-1] / W_top[t, 1:-1, 1:-1])
            - Fa_W_top_WE[t - 1, 1:-1, 1:-1]
            * (Sa_track_top[t, 1:-1, 1:-1] / W_top[t, 1:-1, 1:-1])
            + Fa_W_top_EW[t - 1, 1:-1, 1:-1]
            * (Sa_track_W_top[0, 1:-1, 1:-1] / Sa_W_top[0, 1:-1, 1:-1])
            + Fa_N_top_SN[t - 1, 1:-1, 1:-1]
            * (Sa_track_N_top[0, 1:-1, 1:-1] / Sa_N_top[0, 1:-1, 1:-1])
            - Fa_N_top_NS[t - 1, 1:-1, 1:-1]
            * (Sa_track_top[t, 1:-1, 1:-1] / W_top[t, 1:-1, 1:-1])
            - Fa_S_top_SN[t - 1, 1:-1, 1:-1]
            * (Sa_track_top[t, 1:-1, 1:-1] / W_top[t, 1:-1, 1:-1])
            + Fa_S_top_NS[t - 1, 1:-1, 1:-1]
            * (Sa_track_S_top[0, 1:-1, 1:-1] / Sa_S_top[0, 1:-1, 1:-1])
            + Fa_downward[t - 1, 1:-1, 1:-1]
            * (Sa_track_down[t, 1:-1, 1:-1] / W_down[t, 1:-1, 1:-1])
            - Fa_upward[t - 1, 1:-1, 1:-1]
            * (Sa_track_top[t, 1:-1, 1:-1] / W_top[t, 1:-1, 1:-1])
        )

        # losses to the north and south
        north_loss[t - 1, 0, :] = Fa_N_top_NS[t - 1, 1, :] * (
            Sa_track_top[t, 1, :] / W_top[t, 1, :]
        ) + Fa_N_down_NS[t - 1, 1, :] * (Sa_track_down[t, 1, :] / W_down[t, 1, :])
        south_loss[t - 1, 0, :] = Fa_S_top_SN[t - 1, -2, :] * (
            Sa_track_top[t, -2, :] / W_top[t, -2, :]
        ) + Fa_S_down_SN[t - 1, -2, :] * (Sa_track_down[t, -2, :] / W_down[t, -2, :])

        # Imme added: losses to the east and west
        east_loss[t - 1, 0, :] = Fa_E_top_EW[t - 1, :, -2] * (
            Sa_track_top[t, :, -2] / W_top[t, :, -2]
        ) + Fa_E_down_EW[t - 1, :, -2] * (Sa_track_down[t, :, -2] / W_down[t, :, -2])
        west_loss[t - 1, 0, :] = Fa_W_top_WE[t - 1, :, 1] * (
            Sa_track_top[t, :, 1] / W_top[t, :, 1]
        ) + Fa_W_down_WE[t - 1, :, 1] * (Sa_track_down[t, :, 1] / W_down[t, :, 1])

        # down: add precipitation and subtract evaporation
        Sa_track_after_Fa_P_E_down[0, 1:-1, 1:-1] = (
            Sa_track_after_Fa_down[0, 1:-1, 1:-1]
            + P_region[t - 1, 1:-1, 1:-1] * (W_down[t, 1:-1, 1:-1] / W[t, 1:-1, 1:-1])
            - E[t - 1, 1:-1, 1:-1]
            * (Sa_track_down[t, 1:-1, 1:-1] / W_down[t, 1:-1, 1:-1])
        )

        # top: add precipitation
        Sa_track_after_Fa_P_E_top[0, 1:-1, 1:-1] = Sa_track_after_Fa_top[
            0, 1:-1, 1:-1
        ] + P_region[t - 1, 1:-1, 1:-1] * (W_top[t, 1:-1, 1:-1] / W[t, 1:-1, 1:-1])

        # down and top: redistribute unaccounted water that is otherwise lost from the sytem
        down_to_top[t - 1, :, :] = np.reshape(
            np.maximum(
                0,
                np.reshape(
                    Sa_track_after_Fa_P_E_down, (np.size(Sa_track_after_Fa_P_E_down))
                )
                - np.reshape(W_down[t - 1, :, :], (np.size(W_down[t - 1, :, :]))),
            ),
            (len(latitude), len(longitude)),
        )  # Imme: invalid value encountered in maximum
        top_to_down[t - 1, :, :] = np.reshape(
            np.maximum(
                0,
                np.reshape(
                    Sa_track_after_Fa_P_E_top, (np.size(Sa_track_after_Fa_P_E_top))
                )
                - np.reshape(W_top[t - 1, :, :], (np.size(W_top[t - 1, :, :]))),
            ),
            (len(latitude), len(longitude)),
        )  # Imme: invalid value encountered in maximum
        Sa_track_down[t - 1, 1:-1, 1:-1] = (
            Sa_track_after_Fa_P_E_down[0, 1:-1, 1:-1]
            - down_to_top[t - 1, 1:-1, 1:-1]
            + top_to_down[t - 1, 1:-1, 1:-1]
        )
        Sa_track_top[t - 1, 1:-1, 1:-1] = (
            Sa_track_after_Fa_P_E_top[0, 1:-1, 1:-1]
            - top_to_down[t - 1, 1:-1, 1:-1]
            + down_to_top[t - 1, 1:-1, 1:-1]
        )

    return (
        Sa_track_top,
        Sa_track_down,
        north_loss,
        south_loss,
        east_loss,
        west_loss,
        down_to_top,
        top_to_down,
        water_lost,
    )


# Code
# def get_Sa_track_backward_TIME(latitude,longitude,count_time,divt,timestep,Kvf,Region,Fa_E_top,Fa_N_top,Fa_E_down,Fa_N_down,Fa_Vert,E,P,
#                                            W_top,W_down,Sa_track_top_last,Sa_track_down_last,Sa_time_top_last,Sa_time_down_last):

# I didn't use this piece of code in this study

# create empty array for track and time


def create_empty_array(count_time, divt, latitude, longitude, years):
    Sa_time_top = np.zeros(
        (np.int(count_time * divt) + 1, len(latitude), len(longitude))
    )
    Sa_time_down = np.zeros(
        (np.int(count_time * divt) + 1, len(latitude), len(longitude))
    )
    Sa_track_top = np.zeros(
        (np.int(count_time * divt) + 1, len(latitude), len(longitude))
    )
    Sa_track_down = np.zeros(
        (np.int(count_time * divt) + 1, len(latitude), len(longitude))
    )

    np.savez_compressed(
        datapathea[0],
        Sa_track_top_last=Sa_track_top[0, :, :],
        Sa_track_down_last=Sa_track_down[0, :, :],
    )
    np.savez_compressed(
        datapathea[1], Sa_time_top=Sa_time_top, Sa_time_down=Sa_time_down
    )

    #    sio.savemat(datapathea[0], {'Sa_track_top':Sa_track_top,'Sa_track_down':Sa_track_down},do_compression=True)
    #    sio.savemat(datapathea[1], {'Sa_time_top':Sa_time_top,'Sa_time_down':Sa_time_down},do_compression=True)
    return


# Runtime & Results

start1 = timer()

# The two lines below create empty arrays for first runs/initial values are zero.
previous_data_to_load = datelist[1:][0] + dt.timedelta(days=1)
datapathea = data_path_ea(
    previous_data_to_load.year, previous_data_to_load.month, previous_data_to_load.day
)  # define paths for empty arrays
if veryfirstrun == 1:
    create_empty_array(
        count_time, divt, latitude, longitude, years
    )  # creates empty arrays for first day run
    # so in this specific case for 2011 0 an empty array is created with zeros

for date in datelist[1:]:
    start = timer()

    a = date.day
    yearnumber = date.year
    monthnumber = date.month

    previous_data_to_load = date + dt.timedelta(days=1)
    datapath = data_path(previous_data_to_load, yearnumber, monthnumber, a)

    print(date, previous_data_to_load)
    print(datapath[0])
    # Imme: Hier laad je de getrackte data van de laatste tijdstap, als de laatste tijdstap er neit was dan is die aangemaakt met create_empty_array en zit ie vol met zeros
    loading_ST = np.load(datapath[0])
    # Sa_track_top = loading_ST['Sa_track_top'] # array with zeros #Imme moeten dit zeros zijn of al ingevulde data
    # Sa_track_down = loading_ST['Sa_track_down']
    Sa_track_top_last_1 = loading_ST["Sa_track_top_last"]  # Sa_track_top[0,:,:]
    Sa_track_down_last_1 = loading_ST["Sa_track_down_last"]  # Sa_track_down[0,:,:]
    Sa_track_top_last = np.reshape(
        Sa_track_top_last_1, (1, len(latitude), len(longitude))
    )  # in deze array staan nan en volgens mij hoort dat niet!!
    Sa_track_down_last = np.reshape(
        Sa_track_down_last_1, (1, len(latitude), len(longitude))
    )  # in deze array staan nan en volgens mij hoort dat niet!!

    loading_FS = sio.loadmat(datapath[1], verify_compressed_data_integrity=False)
    Fa_E_top = loading_FS["Fa_E_top"]
    Fa_N_top = loading_FS["Fa_N_top"]
    Fa_E_down = loading_FS["Fa_E_down"]
    Fa_N_down = loading_FS["Fa_N_down"]
    E = loading_FS["E"]
    P = loading_FS["P"]
    W_top = loading_FS["W_top"]
    W_down = loading_FS["W_down"]
    Fa_Vert = loading_FS["Fa_Vert"]

    # call the backward tracking function
    if timetracking == 0:  # I use timetracking = 0
        (
            Sa_track_top,
            Sa_track_down,
            north_loss,
            south_loss,
            east_loss,
            west_loss,
            down_to_top,
            top_to_down,
            water_lost,
        ) = get_Sa_track_backward(
            latitude,
            longitude,
            count_time,
            divt,
            Kvf,
            Region,
            Fa_E_top,
            Fa_N_top,
            Fa_E_down,
            Fa_N_down,
            Fa_Vert,
            E,
            P,
            W_top,
            W_down,
            Sa_track_top_last,
            Sa_track_down_last,
        )
    #    elif timetracking == 1:
    #        loading_STT = sio.loadmat(datapath[2],verify_compressed_data_integrity=False)
    #        Sa_time_top = loading_STT['Sa_time_top'] # [seconds]
    #        Sa_time_down = loading_STT['Sa_time_down']
    #        Sa_time_top_last_1 = Sa_time_top[0,:,:]
    #        Sa_time_down_last_1 = Sa_time_down[0,:,:]
    #        Sa_time_top_last =  np.reshape(Sa_time_top_last_1, (1,len(latitude),len(longitude)))
    #        Sa_time_down_last =  np.reshape(Sa_time_down_last_1, (1,len(latitude),len(longitude)))
    #
    #        Sa_time_top,Sa_time_down,Sa_track_top,Sa_track_down,north_loss,south_loss,down_to_top,top_to_down,water_lost = get_Sa_track_backward_TIME(latitude,longitude,count_time,divt,
    #                                        timestep,Kvf,Region,Fa_E_top,Fa_N_top,Fa_E_down,Fa_N_down,Fa_Vert,E,P,W_top,W_down,Sa_track_top_last,Sa_track_down_last,Sa_time_top_last,Sa_time_down_last)
    # save this data
    # sio.savemat(datapath[3], {'Sa_track_top':Sa_track_top,'Sa_track_down':Sa_track_down,'north_loss':north_loss, 'south_loss':south_loss,'down_to_top':down_to_top,'top_to_down':top_to_down,'water_lost':water_lost},do_compression=True)

    # compute tracked evaporation
    E_track = E[:, :, :] * (Sa_track_down[1:, :, :] / W_down[1:, :, :])

    # save per day
    E_per_day = np.sum(E, axis=0)
    E_track_per_day = np.sum(E_track, axis=0)
    P_per_day = np.sum(P, axis=0)
    Sa_track_down_per_day = np.mean(Sa_track_down[1:, :, :], axis=0)
    Sa_track_top_per_day = np.mean(Sa_track_top[1:, :, :], axis=0)
    W_down_per_day = np.mean(W_down[1:, :, :], axis=0)
    W_top_per_day = np.mean(W_top[1:, :, :], axis=0)

    north_loss_per_day = np.sum(north_loss, axis=0)
    south_loss_per_day = np.sum(south_loss, axis=0)
    east_loss_per_day = np.sum(east_loss, axis=0)
    west_loss_per_day = np.sum(west_loss, axis=0)
    down_to_top_per_day = np.sum(down_to_top, axis=0)
    top_to_down_per_day = np.sum(top_to_down, axis=0)
    water_lost_per_day = np.sum(water_lost, axis=0)

    np.savez_compressed(
        datapath[3],
        Sa_track_top_last=Sa_track_top[0, :, :],
        Sa_track_down_last=Sa_track_down[0, :, :],
        E_per_day=E_per_day,
        E_track_per_day=E_track_per_day,
        P_per_day=P_per_day,
        Sa_track_top_per_day=Sa_track_top_per_day,
        Sa_track_down_per_day=Sa_track_down_per_day,
        W_down_per_day=W_down_per_day,
        W_top_per_day=W_top_per_day,
        north_loss_per_day=north_loss_per_day,
        south_loss_per_day=south_loss_per_day,
        east_loss_per_day=east_loss_per_day,
        west_loss_per_day=west_loss_per_day,
        water_lost_per_day=water_lost_per_day,
    )

    #    if timetracking == 1:
    #        sio.savemat(datapath[4], {'Sa_time_top':Sa_time_top,'Sa_time_down':Sa_time_down},do_compression=True)

    end = timer()
    print(
        "Runtime Sa_track for day "
        + str(a)
        + " in month "
        + str(monthnumber)
        + " in year "
        + str(yearnumber)
        + " is",
        (end - start),
        " seconds.",
    )

end1 = timer()
print("The total runtime of Backtrack_Masterscript is", (end1 - start1), " seconds.")
