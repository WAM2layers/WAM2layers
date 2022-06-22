import datetime as dt
from re import X
import yaml

import numpy as np
from pathlib import Path
import pandas as pd

# Read case configuration
with open("cases/era5_2013.yaml") as f:
    config = yaml.safe_load(f)


datelist = pd.date_range(
    start=config["track_start_date"], end=config["track_end_date"], freq="d", inclusive="left"
)

input_dir = Path(config['preprocessed_data_folder'])
output_dir = Path(config['output_folder']) / "backtrack"

# Check if input dir exists
if not input_dir.exists():
    raise ValueError("Please create the preprocessed_data_folder before running the script")

# Create output dir if it doesn't exist yet
if not output_dir.exists():
    output_dir.mkdir(parents=True)


def input_path(date):
    return f"{input_dir}_{date.strftime('%Y-%m-%d')}_fluxes_storages.nc"


def output_path(date):
    return f"{output_dir}_{date.strftime('%Y-%m-%d')}_Sa_track.nc"


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
        config["interdata_folder"],
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
    Fa_E_upper,
    Fa_N_upper,
    Fa_E_lower,
    Fa_N_lower,
    Fa_Vert,
    E,
    P,
    W_upper,
    W_lower,
    Sa_track_upper_last,
    Sa_track_lower_last,
):

    # make P_region matrix
    Region3D = np.tile(
        np.reshape(Region, [1, len(latitude), len(longitude)]), [len(P[:, 0, 0]), 1, 1]
    )
    P_region = Region3D * P

    # Total moisture in the column
    W = W_upper + W_lower

    # separate the direction of the vertical flux and make it absolute
    Fa_upward = np.zeros(np.shape(Fa_Vert))
    Fa_upward[Fa_Vert <= 0] = Fa_Vert[
        Fa_Vert <= 0
    ]  # in 4th timestep: __main__:1: RuntimeWarning: invalid value encountered in less_equal # not a problem
    Fa_lowerward = np.zeros(np.shape(Fa_Vert))
    Fa_lowerward[Fa_Vert >= 0] = Fa_Vert[Fa_Vert >= 0]
    Fa_upward = np.abs(Fa_upward)

    # include the vertical dispersion
    if Kvf == 0:
        pass
        # do nothing
    else:
        Fa_upward = (1.0 + Kvf) * Fa_upward
        Fa_upward[Fa_Vert >= 0] = Fa_Vert[Fa_Vert >= 0] * Kvf
        Fa_lowerward = (1.0 + Kvf) * Fa_lowerward
        Fa_lowerward[Fa_Vert <= 0] = np.abs(Fa_Vert[Fa_Vert <= 0]) * Kvf

    # define the horizontal fluxes over the boundaries
    # fluxes over the eastern boundary
    Fa_E_upper_boundary = np.zeros(np.shape(Fa_E_upper))
    Fa_E_upper_boundary[:, :, :-1] = 0.5 * (Fa_E_upper[:, :, :-1] + Fa_E_upper[:, :, 1:])
    if config["isglobal"]:
        Fa_E_upper_boundary[:, :, -1] = 0.5 * (Fa_E_upper[:, :, -1] + Fa_E_upper[:, :, 0])
    Fa_E_lower_boundary = np.zeros(np.shape(Fa_E_lower))
    Fa_E_lower_boundary[:, :, :-1] = 0.5 * (Fa_E_lower[:, :, :-1] + Fa_E_lower[:, :, 1:])
    if config["isglobal"]:
        Fa_E_lower_boundary[:, :, -1] = 0.5 * (Fa_E_lower[:, :, -1] + Fa_E_lower[:, :, 0])

    # find out where the positive and negative fluxes are
    Fa_E_upper_pos = np.ones(np.shape(Fa_E_upper))
    Fa_E_lower_pos = np.ones(np.shape(Fa_E_lower))
    Fa_E_upper_pos[Fa_E_upper_boundary < 0] = 0
    Fa_E_lower_pos[Fa_E_lower_boundary < 0] = 0
    Fa_E_upper_neg = Fa_E_upper_pos - 1
    Fa_E_lower_neg = Fa_E_lower_pos - 1

    # separate directions west-east (all positive numbers)
    Fa_E_upper_WE = Fa_E_upper_boundary * Fa_E_upper_pos
    Fa_E_upper_EW = Fa_E_upper_boundary * Fa_E_upper_neg
    Fa_E_lower_WE = Fa_E_lower_boundary * Fa_E_lower_pos
    Fa_E_lower_EW = Fa_E_lower_boundary * Fa_E_lower_neg

    # fluxes over the western boundary
    Fa_W_upper_WE = np.nan * np.zeros(np.shape(P))
    Fa_W_upper_WE[:, :, 1:] = Fa_E_upper_WE[:, :, :-1]
    Fa_W_upper_WE[:, :, 0] = Fa_E_upper_WE[:, :, -1]
    Fa_W_upper_EW = np.nan * np.zeros(np.shape(P))
    Fa_W_upper_EW[:, :, 1:] = Fa_E_upper_EW[:, :, :-1]
    Fa_W_upper_EW[:, :, 0] = Fa_E_upper_EW[:, :, -1]
    Fa_W_lower_WE = np.nan * np.zeros(np.shape(P))
    Fa_W_lower_WE[:, :, 1:] = Fa_E_lower_WE[:, :, :-1]
    Fa_W_lower_WE[:, :, 0] = Fa_E_lower_WE[:, :, -1]
    Fa_W_lower_EW = np.nan * np.zeros(np.shape(P))
    Fa_W_lower_EW[:, :, 1:] = Fa_E_lower_EW[:, :, :-1]
    Fa_W_lower_EW[:, :, 0] = Fa_E_lower_EW[:, :, -1]

    # fluxes over the northern boundary
    # Fa_N_upper_boundary = np.nan*np.zeros(np.shape(Fa_N_upper)); #Imme: why do you multiply here your zeros with nans ?!?!?!? you do not do that for Fa_E_upper_boundary
    # Adapted by Imme
    Fa_N_upper_boundary = np.zeros(np.shape(Fa_N_upper))
    # Imme: why do you multiply here your zeros with nans ?!?!?!? you do not do that for Fa_E_upper_boundary
    Fa_N_upper_boundary[:, 1:, :] = 0.5 * (Fa_N_upper[:, :-1, :] + Fa_N_upper[:, 1:, :])
    # Fa_N_lower_boundary = np.nan*np.zeros(np.shape(Fa_N_lower));
    # Adapted by Imme
    Fa_N_lower_boundary = np.zeros(np.shape(Fa_N_lower))
    # als je niet met np.nan multiplied krijg je in de volgende alinea geen invalid value encountered in less
    # verandert er verder nog wat!??!
    Fa_N_lower_boundary[:, 1:, :] = 0.5 * (Fa_N_lower[:, :-1, :] + Fa_N_lower[:, 1:, :])

    # find out where the positive and negative fluxes are
    Fa_N_upper_pos = np.ones(np.shape(Fa_N_upper))
    Fa_N_lower_pos = np.ones(np.shape(Fa_N_lower))
    Fa_N_upper_pos[
        Fa_N_upper_boundary < 0
    ] = 0  # Invalid value encountered in less, omdat er nan in Fa_N_upper_boundary staan
    Fa_N_lower_pos[
        Fa_N_lower_boundary < 0
    ] = 0  # Invalid value encountered in less, omdat er nan in Fa_N_upper_boundary staan
    Fa_N_upper_neg = Fa_N_upper_pos - 1
    Fa_N_lower_neg = Fa_N_lower_pos - 1

    # separate directions south-north (all positive numbers)
    Fa_N_upper_SN = Fa_N_upper_boundary * Fa_N_upper_pos
    Fa_N_upper_NS = Fa_N_upper_boundary * Fa_N_upper_neg
    Fa_N_lower_SN = Fa_N_lower_boundary * Fa_N_lower_pos
    Fa_N_lower_NS = Fa_N_lower_boundary * Fa_N_lower_neg

    # fluxes over the southern boundary
    Fa_S_upper_SN = np.nan * np.zeros(np.shape(P))
    Fa_S_upper_SN[:, :-1, :] = Fa_N_upper_SN[:, 1:, :]
    Fa_S_upper_NS = np.nan * np.zeros(np.shape(P))
    Fa_S_upper_NS[:, :-1, :] = Fa_N_upper_NS[:, 1:, :]
    Fa_S_lower_SN = np.nan * np.zeros(np.shape(P))
    Fa_S_lower_SN[:, :-1, :] = Fa_N_lower_SN[:, 1:, :]
    Fa_S_lower_NS = np.nan * np.zeros(np.shape(P))
    Fa_S_lower_NS[:, :-1, :] = Fa_N_lower_NS[:, 1:, :]

    # defining size of output
    Sa_track_lower = np.zeros(np.shape(W_lower))
    Sa_track_upper = np.zeros(np.shape(W_upper))

    # assign begin values of output == last (but first index) values of the previous time slot
    Sa_track_lower[-1, :, :] = Sa_track_lower_last
    Sa_track_upper[-1, :, :] = Sa_track_upper_last

    # defining sizes of tracked moisture
    Sa_track_after_Fa_lower = np.zeros(np.shape(Sa_track_lower_last))
    Sa_track_after_Fa_P_E_lower = np.zeros(np.shape(Sa_track_lower_last))
    Sa_track_E_lower = np.zeros(np.shape(Sa_track_lower_last))
    Sa_track_W_lower = np.zeros(np.shape(Sa_track_lower_last))
    Sa_track_N_lower = np.zeros(np.shape(Sa_track_lower_last))
    Sa_track_S_lower = np.zeros(np.shape(Sa_track_lower_last))
    Sa_track_after_Fa_upper = np.zeros(np.shape(Sa_track_upper_last))
    Sa_track_after_Fa_P_E_upper = np.zeros(np.shape(Sa_track_upper_last))
    Sa_track_E_upper = np.zeros(np.shape(Sa_track_upper_last))
    Sa_track_W_upper = np.zeros(np.shape(Sa_track_upper_last))
    Sa_track_N_upper = np.zeros(np.shape(Sa_track_upper_last))
    Sa_track_S_upper = np.zeros(np.shape(Sa_track_upper_last))

    # define sizes of total moisture
    Sa_E_lower = np.zeros(np.shape(Sa_track_lower_last))
    Sa_W_lower = np.zeros(np.shape(Sa_track_lower_last))
    Sa_N_lower = np.zeros(np.shape(Sa_track_lower_last))
    Sa_S_lower = np.zeros(np.shape(Sa_track_lower_last))
    Sa_E_upper = np.zeros(np.shape(Sa_track_upper_last))
    Sa_W_upper = np.zeros(np.shape(Sa_track_upper_last))
    Sa_N_upper = np.zeros(np.shape(Sa_track_upper_last))
    Sa_S_upper = np.zeros(np.shape(Sa_track_upper_last))

    # define variables that find out what happens to the water
    north_loss = np.zeros((np.int(count_time * divt), 1, len(longitude)))
    south_loss = np.zeros((np.int(count_time * divt), 1, len(longitude)))
    east_loss = np.zeros((np.int(count_time * divt), 1, len(latitude)))
    west_loss = np.zeros((np.int(count_time * divt), 1, len(latitude)))
    down_to_upper = np.zeros(np.shape(P))
    top_to_lower = np.zeros(np.shape(P))
    water_lost = np.zeros(np.shape(P))
    water_lost_lower = np.zeros(np.shape(P))
    water_lost_upper = np.zeros(np.shape(P))

    # Sa calculation backward in time
    for t in np.arange(np.int(count_time * divt), 0, -1):
        # down: define values of total moisture
        Sa_E_lower[0, :, :-1] = W_lower[
            t, :, 1:
        ]  # Atmospheric storage of the cell to the east [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors
        # Sa_E_lower[0,:,-1] = W_lower[t,:,0] # Atmospheric storage of the cell to the east [m3]
        Sa_W_lower[0, :, 1:] = W_lower[
            t, :, :-1
        ]  # Atmospheric storage of the cell to the west [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors
        # Sa_W_lower[0,:,0] = W_lower[t,:,-1] # Atmospheric storage of the cell to the west [m3]
        Sa_N_lower[0, 1:, :] = W_lower[
            t, 0:-1, :
        ]  # Atmospheric storage of the cell to the north [m3]
        Sa_S_lower[0, :-1, :] = W_lower[
            t, 1:, :
        ]  # Atmospheric storage of the cell to the south [m3]

        # top: define values of total moisture
        Sa_E_upper[0, :, :-1] = W_upper[
            t, :, 1:
        ]  # Atmospheric storage of the cell to the east [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors
        # Sa_E_upper[0,:,-1] = W_upper[t,:,0] # Atmospheric storage of the cell to the east [m3]
        Sa_W_upper[0, :, 1:] = W_upper[
            t, :, :-1
        ]  # Atmospheric storage of the cell to the west [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors
        # Sa_W_upper[0,:,0] = W_upper[t,:,-1] # Atmospheric storage of the cell to the west [m3]
        Sa_N_upper[0, 1:, :] = W_upper[
            t, :-1, :
        ]  # Atmospheric storage of the cell to the north [m3]
        Sa_S_upper[0, :-1, :] = W_upper[
            t, 1:, :
        ]  # Atmospheric storage of the cell to the south [m3]

        # down: define values of tracked moisture of neighbouring grid cells
        Sa_track_E_lower[0, :, :-1] = Sa_track_lower[
            t, :, 1:
        ]  # Atmospheric tracked storage of the cell to the east [m3]
        if config["isglobal"]:
            Sa_track_E_lower[0, :, -1] = Sa_track_lower[
                t, :, 0
            ]  # Atmospheric tracked storage of the cell to the east [m3]
        Sa_track_W_lower[0, :, 1:] = Sa_track_lower[
            t, :, :-1
        ]  # Atmospheric storage of the cell to the west [m3]
        if config["isglobal"]:
            Sa_track_W_lower[0, :, 0] = Sa_track_lower[
                t, :, -1
            ]  # Atmospheric storage of the cell to the west [m3]
        Sa_track_N_lower[0, 1:, :] = Sa_track_lower[
            t, :-1, :
        ]  # Atmospheric storage of the cell to the north [m3]
        Sa_track_S_lower[0, :-1, :] = Sa_track_lower[
            t, 1:, :
        ]  # Atmospheric storage of the cell to the south [m3]

        # down: calculate with moisture fluxes
        Sa_track_after_Fa_lower[0, 1:-1, 1:-1] = (
            Sa_track_lower[t, 1:-1, 1:-1]
            + Fa_E_lower_WE[t - 1, 1:-1, 1:-1]
            * (Sa_track_E_lower[0, 1:-1, 1:-1] / Sa_E_lower[0, 1:-1, 1:-1])
            - Fa_E_lower_EW[t - 1, 1:-1, 1:-1]
            * (Sa_track_lower[t, 1:-1, 1:-1] / W_lower[t, 1:-1, 1:-1])
            - Fa_W_lower_WE[t - 1, 1:-1, 1:-1]
            * (Sa_track_lower[t, 1:-1, 1:-1] / W_lower[t, 1:-1, 1:-1])
            + Fa_W_lower_EW[t - 1, 1:-1, 1:-1]
            * (Sa_track_W_lower[0, 1:-1, 1:-1] / Sa_W_lower[0, 1:-1, 1:-1])
            + Fa_N_lower_SN[t - 1, 1:-1, 1:-1]
            * (Sa_track_N_lower[0, 1:-1, 1:-1] / Sa_N_lower[0, 1:-1, 1:-1])
            - Fa_N_lower_NS[t - 1, 1:-1, 1:-1]
            * (Sa_track_lower[t, 1:-1, 1:-1] / W_lower[t, 1:-1, 1:-1])
            - Fa_S_lower_SN[t - 1, 1:-1, 1:-1]
            * (Sa_track_lower[t, 1:-1, 1:-1] / W_lower[t, 1:-1, 1:-1])
            + Fa_S_lower_NS[t - 1, 1:-1, 1:-1]
            * (Sa_track_S_lower[0, 1:-1, 1:-1] / Sa_S_lower[0, 1:-1, 1:-1])
            - Fa_lowerward[t - 1, 1:-1, 1:-1]
            * (Sa_track_lower[t, 1:-1, 1:-1] / W_lower[t, 1:-1, 1:-1])
            + Fa_upward[t - 1, 1:-1, 1:-1]
            * (Sa_track_upper[t, 1:-1, 1:-1] / W_upper[t, 1:-1, 1:-1])
        )
        # sometimes you get a Runtimewarning: invalid value encountered in less
        # RuntimeWarning: invalid value encountered in divide
        # ik kan me voorstellen dat die error er komt als er een nan in Sa_track_lower zit en dat je daar dan door moet gaan delen..

        # top: define values of tracked moisture of neighbouring grid cells
        Sa_track_E_upper[0, :, :-1] = Sa_track_upper[
            t, :, 1:
        ]  # Atmospheric tracked storage of the cell to the east [m3]
        if config["isglobal"]:
            Sa_track_E_upper[0, :, -1] = Sa_track_upper[
                t, :, 0
            ]  # Atmospheric tracked storage of the cell to the east [m3]
        Sa_track_W_upper[0, :, 1:] = Sa_track_upper[
            t, :, :-1
        ]  # Atmospheric tracked storage of the cell to the west [m3]
        if config["isglobal"]:
            Sa_track_W_upper[0, :, 0] = Sa_track_upper[
                t, :, -1
            ]  # Atmospheric tracked storage of the cell to the west [m3]
        Sa_track_N_upper[0, 1:, :] = Sa_track_upper[
            t, :-1, :
        ]  # Atmospheric tracked storage of the cell to the north [m3]
        Sa_track_S_upper[0, :-1, :] = Sa_track_upper[
            t, 1:, :
        ]  # Atmospheric tracked storage of the cell to the south [m3]

        # top: calculate with moisture fluxes
        Sa_track_after_Fa_upper[0, 1:-1, 1:-1] = (
            Sa_track_upper[t, 1:-1, 1:-1]
            + Fa_E_upper_WE[t - 1, 1:-1, 1:-1]
            * (Sa_track_E_upper[0, 1:-1, 1:-1] / Sa_E_upper[0, 1:-1, 1:-1])
            - Fa_E_upper_EW[t - 1, 1:-1, 1:-1]
            * (Sa_track_upper[t, 1:-1, 1:-1] / W_upper[t, 1:-1, 1:-1])
            - Fa_W_upper_WE[t - 1, 1:-1, 1:-1]
            * (Sa_track_upper[t, 1:-1, 1:-1] / W_upper[t, 1:-1, 1:-1])
            + Fa_W_upper_EW[t - 1, 1:-1, 1:-1]
            * (Sa_track_W_upper[0, 1:-1, 1:-1] / Sa_W_upper[0, 1:-1, 1:-1])
            + Fa_N_upper_SN[t - 1, 1:-1, 1:-1]
            * (Sa_track_N_upper[0, 1:-1, 1:-1] / Sa_N_upper[0, 1:-1, 1:-1])
            - Fa_N_upper_NS[t - 1, 1:-1, 1:-1]
            * (Sa_track_upper[t, 1:-1, 1:-1] / W_upper[t, 1:-1, 1:-1])
            - Fa_S_upper_SN[t - 1, 1:-1, 1:-1]
            * (Sa_track_upper[t, 1:-1, 1:-1] / W_upper[t, 1:-1, 1:-1])
            + Fa_S_upper_NS[t - 1, 1:-1, 1:-1]
            * (Sa_track_S_upper[0, 1:-1, 1:-1] / Sa_S_upper[0, 1:-1, 1:-1])
            + Fa_lowerward[t - 1, 1:-1, 1:-1]
            * (Sa_track_lower[t, 1:-1, 1:-1] / W_lower[t, 1:-1, 1:-1])
            - Fa_upward[t - 1, 1:-1, 1:-1]
            * (Sa_track_upper[t, 1:-1, 1:-1] / W_upper[t, 1:-1, 1:-1])
        )

        # losses to the north and south
        north_loss[t - 1, 0, :] = Fa_N_upper_NS[t - 1, 1, :] * (
            Sa_track_upper[t, 1, :] / W_upper[t, 1, :]
        ) + Fa_N_lower_NS[t - 1, 1, :] * (Sa_track_lower[t, 1, :] / W_lower[t, 1, :])
        south_loss[t - 1, 0, :] = Fa_S_upper_SN[t - 1, -2, :] * (
            Sa_track_upper[t, -2, :] / W_upper[t, -2, :]
        ) + Fa_S_lower_SN[t - 1, -2, :] * (Sa_track_lower[t, -2, :] / W_lower[t, -2, :])

        # Imme added: losses to the east and west
        east_loss[t - 1, 0, :] = Fa_E_upper_EW[t - 1, :, -2] * (
            Sa_track_upper[t, :, -2] / W_upper[t, :, -2]
        ) + Fa_E_lower_EW[t - 1, :, -2] * (Sa_track_lower[t, :, -2] / W_lower[t, :, -2])
        west_loss[t - 1, 0, :] = Fa_W_upper_WE[t - 1, :, 1] * (
            Sa_track_upper[t, :, 1] / W_upper[t, :, 1]
        ) + Fa_W_lower_WE[t - 1, :, 1] * (Sa_track_lower[t, :, 1] / W_lower[t, :, 1])

        # down: add precipitation and subtract evaporation
        Sa_track_after_Fa_P_E_lower[0, 1:-1, 1:-1] = (
            Sa_track_after_Fa_lower[0, 1:-1, 1:-1]
            + P_region[t - 1, 1:-1, 1:-1] * (W_lower[t, 1:-1, 1:-1] / W[t, 1:-1, 1:-1])
            - E[t - 1, 1:-1, 1:-1]
            * (Sa_track_lower[t, 1:-1, 1:-1] / W_lower[t, 1:-1, 1:-1])
        )

        # top: add precipitation
        Sa_track_after_Fa_P_E_upper[0, 1:-1, 1:-1] = Sa_track_after_Fa_upper[
            0, 1:-1, 1:-1
        ] + P_region[t - 1, 1:-1, 1:-1] * (W_upper[t, 1:-1, 1:-1] / W[t, 1:-1, 1:-1])

        # down and top: redistribute unaccounted water that is otherwise lost from the sytem
        down_to_upper[t - 1, :, :] = np.reshape(
            np.maximum(
                0,
                np.reshape(
                    Sa_track_after_Fa_P_E_lower, (np.size(Sa_track_after_Fa_P_E_lower))
                )
                - np.reshape(W_lower[t - 1, :, :], (np.size(W_lower[t - 1, :, :]))),
            ),
            (len(latitude), len(longitude)),
        )  # Imme: invalid value encountered in maximum
        top_to_lower[t - 1, :, :] = np.reshape(
            np.maximum(
                0,
                np.reshape(
                    Sa_track_after_Fa_P_E_upper, (np.size(Sa_track_after_Fa_P_E_upper))
                )
                - np.reshape(W_upper[t - 1, :, :], (np.size(W_upper[t - 1, :, :]))),
            ),
            (len(latitude), len(longitude)),
        )  # Imme: invalid value encountered in maximum
        Sa_track_lower[t - 1, 1:-1, 1:-1] = (
            Sa_track_after_Fa_P_E_lower[0, 1:-1, 1:-1]
            - down_to_upper[t - 1, 1:-1, 1:-1]
            + top_to_lower[t - 1, 1:-1, 1:-1]
        )
        Sa_track_upper[t - 1, 1:-1, 1:-1] = (
            Sa_track_after_Fa_P_E_upper[0, 1:-1, 1:-1]
            - top_to_lower[t - 1, 1:-1, 1:-1]
            + down_to_upper[t - 1, 1:-1, 1:-1]
        )

    return (
        Sa_track_upper,
        Sa_track_lower,
        north_loss,
        south_loss,
        east_loss,
        west_loss,
        down_to_upper,
        top_to_lower,
        water_lost,
    )


def create_empty_array(count_time, divt, latitude, longitude, ):



# Runtime & Results

# The two lines below create empty arrays for first runs/initial values are zero.
previous_data_to_load = datelist[-1] + dt.timedelta(days=1)
datapathea = data_path_ea(previous_data_to_load)  # define paths for empty arrays

# Load 1 dataset to get grid info
ds = xr.open_dataset(input_path(datelist[0]))


if config["veryfirstrun"]:
    # Create initial field and store to disk
    Sa_track_upper = np.zeros_like(ds.w_upper)
    Sa_track_lower = np.zeros_like(ds.w_upper)

    xr.Dataset(
        {
            "Sa_track_upper": (["time_states", "lat", "lon"], Sa_track_upper),
            "Sa_track_lower": (["time_states", "lat", "lon"], Sa_track_lower),
        },
    ).to_netcdf(datapathea)


for date in datelist[1:]:


    a = date.day
    yearnumber = date.year
    monthnumber = date.month

    previous_data_to_load = date + dt.timedelta(days=1)
    datapath = data_path(previous_data_to_load, yearnumber, monthnumber, a)

    print(date, previous_data_to_load)
    print(datapath[0])
    # Imme: Hier laad je de getrackte data van de laatste tijdstap, als de laatste tijdstap er neit was dan is die aangemaakt met create_empty_array en zit ie vol met zeros
    loading_ST = np.load(datapath[0])
    # Sa_track_upper = loading_ST['Sa_track_upper'] # array with zeros #Imme moeten dit zeros zijn of al ingevulde data
    # Sa_track_lower = loading_ST['Sa_track_lower']
    Sa_track_upper_last_1 = loading_ST["Sa_track_upper_last"]  # Sa_track_upper[0,:,:]
    Sa_track_lower_last_1 = loading_ST["Sa_track_lower_last"]  # Sa_track_lower[0,:,:]
    Sa_track_upper_last = np.reshape(
        Sa_track_upper_last_1, (1, len(latitude), len(longitude))
    )  # in deze array staan nan en volgens mij hoort dat niet!!
    Sa_track_lower_last = np.reshape(
        Sa_track_lower_last_1, (1, len(latitude), len(longitude))
    )  # in deze array staan nan en volgens mij hoort dat niet!!

    loading_FS = sio.loadmat(datapath[1], verify_compressed_data_integrity=False)
    Fa_E_upper = loading_FS["Fa_E_upper"]
    Fa_N_upper = loading_FS["Fa_N_upper"]
    Fa_E_lower = loading_FS["Fa_E_lower"]
    Fa_N_lower = loading_FS["Fa_N_lower"]
    E = loading_FS["E"]
    P = loading_FS["P"]
    W_upper = loading_FS["W_upper"]
    W_lower = loading_FS["W_lower"]
    Fa_Vert = loading_FS["Fa_Vert"]

    # call the backward tracking function
    if not config["timetracking"]:  # I use timetracking: false
        (
            Sa_track_upper,
            Sa_track_lower,
            north_loss,
            south_loss,
            east_loss,
            west_loss,
            down_to_upper,
            top_to_lower,
            water_lost,
        ) = get_Sa_track_backward(
            latitude,
            longitude,
            count_time,
            divt,
            Kvf,
            Region,
            Fa_E_upper,
            Fa_N_upper,
            Fa_E_lower,
            Fa_N_lower,
            Fa_Vert,
            E,
            P,
            W_upper,
            W_lower,
            Sa_track_upper_last,
            Sa_track_lower_last,
        )

    # compute tracked evaporation
    E_track = E[:, :, :] * (Sa_track_lower[1:, :, :] / W_lower[1:, :, :])

    # save per day
    E_per_day = np.sum(E, axis=0)
    E_track_per_day = np.sum(E_track, axis=0)
    P_per_day = np.sum(P, axis=0)
    Sa_track_lower_per_day = np.mean(Sa_track_lower[1:, :, :], axis=0)
    Sa_track_upper_per_day = np.mean(Sa_track_upper[1:, :, :], axis=0)
    W_lower_per_day = np.mean(W_lower[1:, :, :], axis=0)
    W_upper_per_day = np.mean(W_upper[1:, :, :], axis=0)

    north_loss_per_day = np.sum(north_loss, axis=0)
    south_loss_per_day = np.sum(south_loss, axis=0)
    east_loss_per_day = np.sum(east_loss, axis=0)
    west_loss_per_day = np.sum(west_loss, axis=0)
    down_to_upper_per_day = np.sum(down_to_upper, axis=0)
    top_to_lower_per_day = np.sum(top_to_lower, axis=0)
    water_lost_per_day = np.sum(water_lost, axis=0)

    np.savez_compressed(
        datapath[3],
        Sa_track_upper_last=Sa_track_upper[0, :, :],
        Sa_track_lower_last=Sa_track_lower[0, :, :],
        E_per_day=E_per_day,
        E_track_per_day=E_track_per_day,
        P_per_day=P_per_day,
        Sa_track_upper_per_day=Sa_track_upper_per_day,
        Sa_track_lower_per_day=Sa_track_lower_per_day,
        W_lower_per_day=W_lower_per_day,
        W_upper_per_day=W_upper_per_day,
        north_loss_per_day=north_loss_per_day,
        south_loss_per_day=south_loss_per_day,
        east_loss_per_day=east_loss_per_day,
        west_loss_per_day=west_loss_per_day,
        water_lost_per_day=water_lost_per_day,
    )
