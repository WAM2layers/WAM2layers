import yaml
import xarray as xr

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
    return f"{input_dir}/{date.strftime('%Y-%m-%d')}_fluxes_storages.nc"


def output_path(date):
    return f"{output_dir}/{date.strftime('%Y-%m-%d')}_Sa_track.nc"


def backtrack(
    latitude,
    longitude,
    Kvf,
    region,
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
    P_region = region * P

    # Total moisture in the column
    W = W_upper + W_lower

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
    Fa_E_upper_boundary = np.zeros(np.shape(Fa_E_upper))
    Fa_E_upper_boundary[:, :, :-1] = 0.5 * (Fa_E_upper[:, :, :-1] + Fa_E_upper[:, :, 1:])
    if config["periodic_boundary"]:
        Fa_E_upper_boundary[:, :, -1] = 0.5 * (Fa_E_upper[:, :, -1] + Fa_E_upper[:, :, 0])
    Fa_E_lower_boundary = np.zeros(np.shape(Fa_E_lower))
    Fa_E_lower_boundary[:, :, :-1] = 0.5 * (Fa_E_lower[:, :, :-1] + Fa_E_lower[:, :, 1:])
    if config["periodic_boundary"]:
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
    # Imme: why do you multiply here your zeros with nans ?!?!?!? you do not do that for Fa_E_upper_boundary
    # Fa_N_upper_boundary = np.nan*np.zeros(np.shape(Fa_N_upper));
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
    Fa_N_upper_pos[Fa_N_upper_boundary < 0] = 0  # Invalid value encountered in less
    Fa_N_lower_pos[Fa_N_lower_boundary < 0] = 0  # Invalid value encountered in less
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
    ntime = w_upper.shape[0]
    north_loss = np.zeros((ntime, 1, len(longitude)))
    south_loss = np.zeros((ntime, 1, len(longitude)))
    east_loss = np.zeros((ntime, 1, len(latitude)))
    west_loss = np.zeros((ntime, 1, len(latitude)))
    down_to_upper = np.zeros(np.shape(P))
    top_to_lower = np.zeros(np.shape(P))
    water_lost = np.zeros(np.shape(P))

    # Sa calculation backward in time
    for t in np.arange(ntime-1, 0, -1):

        # Upper layer: define values of total moisture
        Sa_N_lower[0, 1:, :] = W_lower[t, 0:-1, :]
        Sa_S_lower[0, :-1, :] = W_lower[t, 1:, :]
        Sa_E_lower[0, :, :-1] = W_lower[t, :, 1:]
        Sa_W_lower[0, :, 1:] = W_lower[t, :, :-1]
        # TODO: fix div by zero issue for periodic_boundary
        # Sa_E_lower[0,:,-1] = W_lower[t,:,0]
        # Sa_W_lower[0,:,0] = W_lower[t,:,-1]

        # Lower layer: define values of total moisture
        Sa_N_upper[0, 1:, :] = W_upper[t, :-1, :]
        Sa_S_upper[0, :-1, :] = W_upper[t, 1:, :]
        Sa_E_upper[0, :, :-1] = W_upper[t, :, 1:]
        Sa_W_upper[0, :, 1:] = W_upper[t, :, :-1]
        # TODO: fix div by zero issue for periodic_boundary
        # Sa_E_upper[0,:,-1] = W_upper[t,:,0]
        # Sa_W_upper[0,:,0] = W_upper[t,:,-1]

        # Lower layer: define values of tracked moisture of neighbouring grid cells
        Sa_track_N_lower[0, 1:, :] = Sa_track_lower[t, :-1, :]
        Sa_track_S_lower[0, :-1, :] = Sa_track_lower[t, 1:, :]
        Sa_track_E_lower[0, :, :-1] = Sa_track_lower[t, :, 1:]
        Sa_track_W_lower[0, :, 1:] = Sa_track_lower[t, :, :-1]
        if config["periodic_boundary"]:
            Sa_track_E_lower[0, :, -1] = Sa_track_lower[t, :, 0]
            Sa_track_W_lower[0, :, 0] = Sa_track_lower[t, :, -1]

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
            - Fa_downward[t - 1, 1:-1, 1:-1]
            * (Sa_track_lower[t, 1:-1, 1:-1] / W_lower[t, 1:-1, 1:-1])
            + Fa_upward[t - 1, 1:-1, 1:-1]
            * (Sa_track_upper[t, 1:-1, 1:-1] / W_upper[t, 1:-1, 1:-1])
        )

        # top: define values of tracked moisture of neighbouring grid cells
        Sa_track_N_upper[0, 1:, :] = Sa_track_upper[t, :-1, :]
        Sa_track_S_upper[0, :-1, :] = Sa_track_upper[t, 1:, :]
        Sa_track_E_upper[0, :, :-1] = Sa_track_upper[t, :, 1:]
        Sa_track_W_upper[0, :, 1:] = Sa_track_upper[t, :, :-1]
        if config["periodic_boundary"]:
            Sa_track_E_upper[0, :, -1] = Sa_track_upper[t, :, 0]
            Sa_track_W_upper[0, :, 0] = Sa_track_upper[t, :, -1]


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
            + Fa_downward[t - 1, 1:-1, 1:-1]
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


if config["veryfirstrun"]:
    # Create initial state (tracked moisture == 0)

    # Load 1 dataset to get grid info
    ds = xr.open_dataset(input_path(datelist[0]))

    s_track_upper = np.zeros_like(ds.w_upper)
    s_track_lower = np.zeros_like(ds.w_upper)

    states = xr.Dataset(
        {
            "s_track_upper": (["time_states", "lat", "lon"], s_track_upper),
            "s_track_lower": (["time_states", "lat", "lon"], s_track_lower),
        }
    ).isel(time_states=[0])

region = xr.open_dataset(config['region']).isel(time=0).lsm.rename(latitude='lat', longitude='lon').values

for date in reversed(datelist):
    print(date)

    s_track_upper = states.s_track_upper.values
    s_track_lower = states.s_track_lower.values

    fluxes = xr.open_dataset(input_path(date))
    fa_e_upper = fluxes["fa_e_upper"].values
    fa_n_upper = fluxes["fa_n_upper"].values
    fa_e_lower = fluxes["fa_e_lower"].values
    fa_n_lower = fluxes["fa_n_lower"].values
    evap = fluxes["evap"].values
    precip = fluxes["precip"].values
    w_upper = fluxes["w_upper"].values
    w_lower = fluxes["w_lower"].values
    fa_vert = fluxes["fa_vert"].values

    (s_track_upper, s_track_lower, north_loss, south_loss, east_loss, west_loss,
        down_to_upper, top_to_lower, water_lost) = backtrack(
        ds.lat,
        ds.lon,
        config['kvf'],
        region,
        fa_e_upper,
        fa_n_upper,
        fa_e_lower,
        fa_n_lower,
        fa_vert,
        evap,
        precip,
        w_upper,
        w_lower,
        s_track_upper,
        s_track_lower,
    )

    # compute tracked evaporation
    e_track = evap[:, :, :] * (s_track_lower[1:, :, :] / w_lower[1:, :, :])

    # Save per day
    e_per_day = np.sum(evap, axis=0)
    e_track_per_day = np.sum(e_track, axis=0)
    p_per_day = np.sum(precip, axis=0)
    s_track_lower_per_day = np.mean(s_track_lower[1:, :, :], axis=0)
    s_track_upper_per_day = np.mean(s_track_upper[1:, :, :], axis=0)
    w_lower_per_day = np.mean(w_lower[1:, :, :], axis=0)
    w_upper_per_day = np.mean(w_upper[1:, :, :], axis=0)

    north_loss_per_day = np.sum(north_loss, axis=0)
    south_loss_per_day = np.sum(south_loss, axis=0)
    east_loss_per_day = np.sum(east_loss, axis=0)
    west_loss_per_day = np.sum(west_loss, axis=0)
    down_to_upper_per_day = np.sum(down_to_upper, axis=0)
    top_to_lower_per_day = np.sum(top_to_lower, axis=0)
    water_lost_per_day = np.sum(water_lost, axis=0)

    # Write output to file
    # TODO: add (and cleanup) coordinates and units
    xr.Dataset(
        {
            "s_track_upper_last": (["lat", "lon"], s_track_upper[0, :, :]),  # Keep last state for a restart
            "s_track_lower_last": (["lat", "lon"], s_track_lower[0, :, :]),
            "e_per_day": (["lat", "lon"], e_per_day),
            "e_track_per_day": (["lat", "lon"], e_track_per_day),
            "p_per_day": (["lat", "lon"], p_per_day),
            "s_track_upper_per_day": (["lat", "lon"], s_track_upper_per_day),
            "s_track_lower_per_day": (["lat", "lon"], s_track_lower_per_day),
            "w_lower_per_day": (["lat", "lon"], w_lower_per_day),
            "w_upper_per_day": (["lat", "lon"], w_upper_per_day),
            "north_loss_per_day": (["lat1", "lon"], north_loss_per_day),
            "south_loss_per_day": (["lat1", "lon"], south_loss_per_day),
            "east_loss_per_day": (["lon1", "lat"], east_loss_per_day),
            "west_loss_per_day": (["lon1", "lat",], west_loss_per_day),
            "water_lost_per_day": (["lat", "lon"], water_lost_per_day),
        }
    ).to_netcdf(output_path(date))
