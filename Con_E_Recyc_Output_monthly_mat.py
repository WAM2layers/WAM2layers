# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 13:24:45 2016
@author: Ent00002
"""

"""
Created on Mon Feb 18 15:30:43 2019

@author: bened003
"""

# This script is almost similar as the Con_E_Recyc_Output script from WAM-2layers from Ruud van der Ent

# We have implemented a datelist function so the model can run for multiple years without having problems with leap years

# Import libraries

import calendar
import datetime
import datetime as dt
import os
from datetime import timedelta
from timeit import default_timer as timer

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


# BEGIN OF INPUT (FILL THIS IN)
months_length_leap = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
months_length_nonleap = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
years = np.arange(
    2002, 2007
)  # fill in the years # If I fill in more than one year than I need to set the months to 12

# Manage the extent of your dataset (FILL THIS IN)
# Define the latitude and longitude cell numbers to consider and corresponding lakes that should be considered part of the land
latnrs = np.arange(0, 267)  # minimal domain
lonnrs = np.arange(0, 444)

# obtain the constants
lsm_data_ECEarth_T799 = "landseamask_ECearth_T799.nc"  # insert landsea mask here

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

interdata_folder = (
    "Interdata_ECEarth/PresentMember2_correct"  # Inster Interdata folder here
)
sub_interdata_folder = os.path.join(
    interdata_folder, "Regional_backward_daily"
)  # Insert sub-interdata folder here
output_folder = "Output_ECEarth/PresentMember2_correct/"  # Insert Output data here

daily = 1
timetracking = 0  # 0 for not tracking time and 1 for tracking time
# END OF INPUT

# Datapaths (FILL THIS IN)


def data_path(y, a, month, years, timetracking):
    load_Sa_track = os.path.join(
        sub_interdata_folder,
        str(y) + "-" + str(month).zfill(2) + "-" + str(a).zfill(2) + "Sa_track.npz",
    )

    load_Sa_time = os.path.join(
        sub_interdata_folder,
        str(y) + "-" + str(month).zfill(2) + "-" + str(a).zfill(2) + "Sa_time.npz",
    )

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
        output_folder,
        "E_track_regional_full"
        + str(years[0])
        + "-"
        + str(years[-1])
        + "-timetracking"
        + str(timetracking),
    )

    save_path_daily = os.path.join(
        output_folder,
        "E_track_regional_daily_full" + str(y) + "-timetracking" + str(timetracking),
    )

    return (
        load_Sa_track,
        load_Sa_time,
        load_fluxes_and_storages,
        save_path,
        save_path_daily,
    )


# Runtime & Results

start1 = timer()
startyear = years[0]

E_per_year_per_month = np.zeros((len(years), 12, len(latitude), len(longitude)))
E_track_per_year_per_month = np.zeros((len(years), 12, len(latitude), len(longitude)))
P_per_year_per_month = np.zeros((len(years), 12, len(latitude), len(longitude)))
Sa_track_down_per_year_per_month = np.zeros(
    (len(years), 12, len(latitude), len(longitude))
)
Sa_track_top_per_year_per_month = np.zeros(
    (len(years), 12, len(latitude), len(longitude))
)
W_down_per_year_per_month = np.zeros((len(years), 12, len(latitude), len(longitude)))
W_top_per_year_per_month = np.zeros((len(years), 12, len(latitude), len(longitude)))
north_loss_per_year_per_month = np.zeros((len(years), 12, 1, len(longitude)))
south_loss_per_year_per_month = np.zeros((len(years), 12, 1, len(longitude)))
east_loss_per_year_per_month = np.zeros((len(years), 12, 1, len(latitude)))
west_loss_per_year_per_month = np.zeros((len(years), 12, 1, len(latitude)))
down_to_top_per_year_per_month = np.zeros(
    (len(years), 12, len(latitude), len(longitude))
)
top_to_down_per_year_per_month = np.zeros(
    (len(years), 12, len(latitude), len(longitude))
)
water_lost_per_year_per_month = np.zeros(
    (len(years), 12, len(latitude), len(longitude))
)

for year in years[:]:
    start = timer()

    ly = int(calendar.isleap(year))
    final_time = 364 + ly

    if year != 2006:  # no leap year
        datelist = get_times_daily(dt.date(year, 1, 1), dt.date(year, 12, 31))
    else:  # leap
        datelist = get_times_daily(dt.date(year, 1, 1), dt.date(year, 12, 30))

    E_per_day = np.zeros((365 + ly, len(latitude), len(longitude)))
    E_track_per_day = np.zeros((365 + ly, len(latitude), len(longitude)))
    P_per_day = np.zeros((365 + ly, len(latitude), len(longitude)))
    Sa_track_down_per_day = np.zeros((365 + ly, len(latitude), len(longitude)))
    Sa_track_top_per_day = np.zeros((365 + ly, len(latitude), len(longitude)))
    W_down_per_day = np.zeros((365 + ly, len(latitude), len(longitude)))
    W_top_per_day = np.zeros((365 + ly, len(latitude), len(longitude)))
    north_loss_per_day = np.zeros((365 + ly, 1, len(longitude)))
    south_loss_per_day = np.zeros((365 + ly, 1, len(longitude)))
    east_loss_per_day = np.zeros((365 + ly, 1, len(latitude)))
    west_loss_per_day = np.zeros((365 + ly, 1, len(latitude)))
    down_to_top_per_day = np.zeros((365 + ly, len(latitude), len(longitude)))
    top_to_down_per_day = np.zeros((365 + ly, len(latitude), len(longitude)))
    water_lost_per_day = np.zeros((365 + ly, len(latitude), len(longitude)))
    # water_lost_top_per_day = np.zeros((365+ly,len(latitude),len(longitude)))

    for i, date in enumerate(datelist):

        a = date.day
        yearnumber = date.year
        monthnumber = date.month
        print(i, yearnumber, monthnumber, a)

        datapath = data_path(yearnumber, a, monthnumber, years, timetracking)

        print(datapath[0])

        if i > final_time:  # a = 365 (366th index) and not a leapyear\
            pass
        else:
            # load tracked data
            loading_ST = np.load(
                datapath[0]
            )  # ,verify_compressed_data_integrity=False)

            # load the total moisture data from fluxes and storages
            loading_FS = sio.loadmat(
                datapath[2], verify_compressed_data_integrity=False
            )

            # save per day
            E_per_day[i, :, :] = loading_ST["E_per_day"]
            E_track_per_day[i, :, :] = loading_ST["E_track_per_day"]
            P_per_day[i, :, :] = loading_ST["P_per_day"]
            Sa_track_down_per_day[i, :, :] = loading_ST["Sa_track_down_per_day"]
            Sa_track_top_per_day[i, :, :] = loading_ST["Sa_track_top_per_day"]
            W_down_per_day[i, :, :] = loading_ST["W_down_per_day"]
            W_top_per_day[i, :, :] = loading_ST["W_top_per_day"]

            north_loss_per_day[i, :, :] = loading_ST["north_loss_per_day"]
            south_loss_per_day[i, :, :] = loading_ST["south_loss_per_day"]
            east_loss_per_day[i, :, :] = loading_ST["east_loss_per_day"]
            west_loss_per_day[i, :, :] = loading_ST["west_loss_per_day"]
            # down_to_top_per_day[i,:,:] = np.sum(down_to_top, axis =0)
            # top_to_down_per_day[i,:,:] = np.sum(top_to_down, axis =0)
            water_lost_per_day[i, :, :] = loading_ST["water_lost_per_day"]

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
                " seconds",
            )

    if daily == 1:
        if timetracking == 0:  # create dummy values
            Sa_time_down_per_day = 0
            Sa_time_top_per_day = 0
            E_time_per_day = 0
        # save per day
        np.savez_compressed(
            datapath[4],
            E_per_day=E_per_day,
            E_track_per_day=E_track_per_day,
            P_per_day=P_per_day,
            Sa_track_down_per_day=Sa_track_down_per_day,
            Sa_track_top_per_day=Sa_track_top_per_day,
            Sa_time_down_per_day=Sa_time_down_per_day,
            Sa_time_top_per_day=Sa_time_top_per_day,
            W_down_per_day=W_down_per_day,
            W_top_per_day=W_top_per_day,
            E_time_per_day=E_time_per_day,
            water_lost_per_day=water_lost_per_day,
        )  # , water_lost_top_per_day=water_lost_top_per_day)#},do_compression=True)

    # values per month
    for m in range(12):
        first_day = int(datetime.date(year, m + 1, 1).strftime("%j"))
        last_day = int(
            datetime.date(year, m + 1, calendar.monthrange(year, m + 1)[1]).strftime(
                "%j"
            )
        )
        days = np.arange(first_day, last_day + 1) - 1  # -1 because Python is zero-based

        E_per_year_per_month[year - startyear, m, :, :] = np.squeeze(
            np.sum(E_per_day[days, :, :], axis=0)
        )
        E_track_per_year_per_month[year - startyear, m, :, :] = np.squeeze(
            np.sum(E_track_per_day[days, :, :], axis=0)
        )
        P_per_year_per_month[year - startyear, m, :, :] = np.squeeze(
            np.sum(P_per_day[days, :, :], axis=0)
        )
        Sa_track_down_per_year_per_month[year - startyear, m, :, :] = np.squeeze(
            np.mean(Sa_track_down_per_day[days, :, :], axis=0)
        )
        Sa_track_top_per_year_per_month[year - startyear, m, :, :] = np.squeeze(
            np.mean(Sa_track_top_per_day[days, :, :], axis=0)
        )
        W_down_per_year_per_month[year - startyear, m, :, :] = np.squeeze(
            np.mean(W_down_per_day[days, :, :], axis=0)
        )
        W_top_per_year_per_month[year - startyear, m, :, :] = np.squeeze(
            np.mean(W_top_per_day[days, :, :], axis=0)
        )
        north_loss_per_year_per_month[year - startyear, m, :, :] = np.squeeze(
            np.sum(north_loss_per_day[days, :, :], axis=0)
        )
        south_loss_per_year_per_month[year - startyear, m, :, :] = np.squeeze(
            np.sum(south_loss_per_day[days, :, :], axis=0)
        )
        east_loss_per_year_per_month[year - startyear, m, :, :] = np.squeeze(
            np.sum(east_loss_per_day[days, :, :], axis=0)
        )
        west_loss_per_year_per_month[year - startyear, m, :, :] = np.squeeze(
            np.sum(west_loss_per_day[days, :, :], axis=0)
        )
        # down_to_top_per_year_per_month[year-startyear,m,:,:] = (np.squeeze(np.sum(down_to_top_per_day[days,:,:], axis = 0)))
        # top_to_down_per_year_per_month[year-startyear,m,:,:] = (np.squeeze(np.sum(top_to_down_per_day[days,:,:], axis = 0)))
        water_lost_per_year_per_month[year - startyear, m, :, :] = np.squeeze(
            np.sum(water_lost_per_day[days, :, :], axis=0)
        )

        # hallo
        if timetracking == 0:
            Sa_time_down_per_year_per_month = 0
            Sa_time_top_per_year_per_month = 0
            E_time_per_year_per_month = 0

# save monthly data
np.savez_compressed(
    datapath[3],
    E_per_year_per_month=E_per_year_per_month,
    E_track_per_year_per_month=E_track_per_year_per_month,
    P_per_year_per_month=P_per_year_per_month,
    Sa_track_down_per_year_per_month=Sa_track_down_per_year_per_month,
    Sa_track_top_per_year_per_month=Sa_track_top_per_year_per_month,
    Sa_time_down_per_year_per_month=Sa_time_down_per_year_per_month,
    Sa_time_top_per_year_per_month=Sa_time_top_per_year_per_month,
    E_time_per_year_per_month=E_time_per_year_per_month,
    W_down_per_year_per_month=W_down_per_year_per_month,
    W_top_per_year_per_month=W_top_per_year_per_month,
    north_loss_per_year_per_month=north_loss_per_year_per_month,
    south_loss_per_year_per_month=south_loss_per_year_per_month,
    east_loss_per_year_per_month=east_loss_per_year_per_month,
    west_loss_per_year_per_month=west_loss_per_year_per_month,
    down_to_top_per_year_per_month=down_to_top_per_year_per_month,
    top_to_down_per_year_per_month=top_to_down_per_year_per_month,
    water_lost_per_year_per_month=water_lost_per_year_per_month,
)  # , water_lost_per_year_per_month=water_lost_top_per_year_per_month)

end1 = timer()
print("The total runtime of Con_E_Recyc_Output is", (end1 - start1), " seconds.")
