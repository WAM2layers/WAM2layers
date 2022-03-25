# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 13:24:45 2016
@author: Ent00002
"""

"""
Created on Mon Feb 18 15:30:43 2019

@author: bened003
"""

# This script is almost similar as the Hor_Fluxes_Output script from WAM-2layers from Ruud van der Ent

# We have implemented a datelist function so the model can run for multiple years without having problems with leap years


# Import libraries

import numpy as np
import scipy.io as sio
import calendar
import datetime
from getconstants_pressure_ECEarth import getconstants_pressure_ECEarth
from timeit import default_timer as timer
import os
from datetime import timedelta
import datetime as dt

# to create datelist
def get_times_daily(startdate, enddate):
    """generate a dictionary with date/times"""
    numdays = enddate - startdate
    dateList = []
    for x in range(0, numdays.days + 1):
        dateList.append(startdate + dt.timedelta(days=x))
    return dateList


# BEGIN OF INPUT (FILL THIS IN)
months_length_leap = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
months_length_nonleap = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
years = np.arange(
    2002, 2007
)  # fill in the years # If I fill in more than one year than I need to set the months to 12
# create datelist

# Manage the extent of your dataset (FILL THIS IN)
# Define the latitude and longitude cell numbers to consider and corresponding lakes that should be considered part of the land
latnrs = np.arange(0, 267)  # minimal domain
lonnrs = np.arange(0, 444)

isglobal = 0  # fill in 1 for global computations (i.e. Earth round), fill in 0 for a local domain with boundaries

daily = 1  # 1 for writing out daily data, 0 for only monthly data

# END OF INPUT
# Datapaths (FILL THIS IN)

lsm_data_ECEarth_T799 = "landseamask_ECearth_T799.nc"  # insert landseamask

interdata_folder = (
    "Interdata_ECEarth/PresentMember4_correct/"  # insert interdata folder here
)
input_folder = "Inputdata_ECEarth/PresentT799Member4/"  # insert input data folder here
output_folder = "Output_ECEarth/PresentMember4_correct"  # insert output folder here


def data_path(y, month, a, years):
    load_fluxes_and_storages = os.path.join(
        interdata_folder,
        str(y)
        + "-"
        + str(month).zfill(2)
        + "-"
        + str(a).zfill(2)
        + "fluxes_storages.mat",
    )

    save_path = os.path.join(
        output_folder, "Hor_Fluxes_full-" + str(years[0]) + "-" + str(years[-1])
    )
    save_path_daily = os.path.join(output_folder, "Hor_Fluxes_daily_full-" + str(y))

    return load_fluxes_and_storages, save_path, save_path_daily


# Runtime % Results

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
) = getconstants_pressure_ECEarth(latnrs, lonnrs, lsm_data_ECEarth_T799)

startyear = years[0]
Fa_E_down_per_year_per_month = np.zeros((len(years), 12, len(latitude), len(longitude)))
Fa_E_top_per_year_per_month = np.zeros((len(years), 12, len(latitude), len(longitude)))
Fa_N_down_per_year_per_month = np.zeros((len(years), 12, len(latitude), len(longitude)))
Fa_N_top_per_year_per_month = np.zeros((len(years), 12, len(latitude), len(longitude)))
Fa_Vert_per_year_per_month = np.zeros((len(years), 12, len(latitude), len(longitude)))

# from other script
for year in years[:]:
    start = timer()

    if year != 2006:  # if no leap year # specific for my dataset as 2006 is a leap year
        datelist = get_times_daily(dt.date(year, 1, 1), dt.date(year, 12, 31))
    else:  # year == 2006 # if leap year
        datelist = get_times_daily(dt.date(year, 1, 1), dt.date(year, 12, 30))

    ly = int(calendar.isleap(year))
    final_time = 364 + ly

    Fa_E_down_per_day = np.zeros((365 + ly, len(latitude), len(longitude)))
    Fa_E_top_per_day = np.zeros((365 + ly, len(latitude), len(longitude)))
    Fa_N_down_per_day = np.zeros((365 + ly, len(latitude), len(longitude)))
    Fa_N_top_per_day = np.zeros((365 + ly, len(latitude), len(longitude)))
    Fa_Vert_per_day = np.zeros((365 + ly, len(latitude), len(longitude)))

    for i, date in enumerate(datelist[:]):

        a = date.day
        yearnumber = date.year
        monthnumber = date.month
        print(i, yearnumber, monthnumber, a)

        datapath = data_path(yearnumber, monthnumber, a, years)

        if i > final_time:  # a = 365 and not a leapyear
            pass
        else:

            # load horizontal fluxes
            loading_FS = sio.loadmat(
                datapath[0], verify_compressed_data_integrity=False
            )
            Fa_E_top = loading_FS["Fa_E_top"]
            Fa_N_top = loading_FS["Fa_N_top"]
            Fa_E_down = loading_FS["Fa_E_down"]
            Fa_N_down = loading_FS["Fa_N_down"]
            Fa_Vert = loading_FS["Fa_Vert"]

            # save per day
            Fa_E_down_per_day[i, :, :] = np.sum(Fa_E_down, axis=0)
            Fa_E_top_per_day[i, :, :] = np.sum(Fa_E_top, axis=0)
            Fa_N_down_per_day[i, :, :] = np.sum(Fa_N_down, axis=0)
            Fa_N_top_per_day[i, :, :] = np.sum(Fa_N_top, axis=0)
            Fa_Vert_per_day[i, :, :] = np.sum(Fa_Vert, axis=0)

            # timer
            end = timer()
            print(
                "Runtime output for day "
                + str(a)
                + "in month "
                + str(monthnumber)
                + " in year "
                + str(yearnumber)
                + " is",
                (end - start),
                " seconds.",
            )

    # save daily fluxes on disk
    if daily == 1:
        np.savez_compressed(
            datapath[2],
            Fa_E_down_per_day=Fa_E_down_per_day,
            Fa_E_top_per_day=Fa_E_top_per_day,
            Fa_N_down_per_day=Fa_N_down_per_day,
            Fa_N_top_per_day=Fa_N_top_per_day,
            Fa_Vert_per_day=Fa_Vert_per_day,
        )  #

    for m in range(12):
        first_day = int(datetime.date(year, m + 1, 1).strftime("%j"))
        last_day = int(
            datetime.date(year, m + 1, calendar.monthrange(year, m + 1)[1]).strftime(
                "%j"
            )
        )
        days = np.arange(first_day, last_day + 1) - 1  # -1 because Python is zero-based
        print(year, first_day, last_day, days, startyear)

        Fa_E_down_per_year_per_month[year - startyear, m, :, :] = np.squeeze(
            np.sum(Fa_E_down_per_day[days, :, :], axis=0)
        )
        Fa_E_top_per_year_per_month[year - startyear, m, :, :] = np.squeeze(
            np.sum(Fa_E_top_per_day[days, :, :], axis=0)
        )
        Fa_N_down_per_year_per_month[year - startyear, m, :, :] = np.squeeze(
            np.sum(Fa_N_down_per_day[days, :, :], axis=0)
        )
        Fa_N_top_per_year_per_month[year - startyear, m, :, :] = np.squeeze(
            np.sum(Fa_N_top_per_day[days, :, :], axis=0)
        )
        Fa_Vert_per_year_per_month[year - startyear, m, :, :] = np.squeeze(
            np.sum(Fa_Vert_per_day[days, :, :], axis=0)
        )

# save monthly fluxes on disk
np.savez_compressed(
    datapath[1],
    Fa_E_down_per_year_per_month=Fa_E_down_per_year_per_month,
    Fa_E_top_per_year_per_month=Fa_E_top_per_year_per_month,
    Fa_N_down_per_year_per_month=Fa_N_down_per_year_per_month,
    Fa_N_top_per_year_per_month=Fa_N_top_per_year_per_month,
    Fa_Vert_per_year_per_month=Fa_Vert_per_year_per_month,
)

end1 = timer()
print("The total runtime is", (end1 - start1), " seconds.")
