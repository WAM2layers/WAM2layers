import yaml
import xarray as xr

import numpy as np
from pathlib import Path
import pandas as pd

# Read case configuration
with open("cases/era5_2013_local.yaml") as f:
    config = yaml.safe_load(f)


datelist = pd.date_range(
    start=config["track_start_date"], end=config["track_end_date"], freq="d", inclusive="left"
)

input_dir = Path(config['preprocessed_data_folder']).expanduser()
output_dir = Path(config['output_folder']).expanduser() / "backtrack"

# Check if input dir exists
if not input_dir.exists():
    raise ValueError("Please create the preprocessed_data_folder before running the script")

# Create output dir if it doesn't exist yet
if not output_dir.exists():
    output_dir.mkdir(parents=True)


def input_path(date):
    return f"{input_dir}/{date.strftime('%Y-%m-%d')}_fluxes_storages.nc"


def output_path(date):
    return f"{output_dir}/{date.strftime('%Y-%m-%d')}_s_track.nc"


def to_boundaries_zonal(fa_e, p):
    """Define the horizontal fluxes over the east/west boundaries."""
    fa_e_boundary = np.zeros_like(fa_e)
    fa_e_boundary[:, :, :-1] = 0.5 * (fa_e[:, :, :-1] + fa_e[:, :, 1:])
    if config["periodic_boundary"]:
        fa_e_boundary[:, :, -1] = 0.5 * (fa_e[:, :, -1] + fa_e[:, :, 0])

    # find out where the positive and negative fluxes are
    fa_e_pos = np.ones_like(fa_e)
    fa_e_neg = fa_e_pos - 1
    fa_e_pos[fa_e_boundary < 0] = 0

    # separate directions west-east (all positive numbers)
    fa_e_we = fa_e_boundary * fa_e_pos
    fa_e_ew = fa_e_boundary * fa_e_neg

    # fluxes over the western boundary
    fa_w_we = np.nan * np.zeros_like(p)
    fa_w_ew = np.nan * np.zeros_like(p)
    fa_w_we[:, :, 1:] = fa_e_we[:, :, :-1]
    fa_w_we[:, :, 0] = fa_e_we[:, :, -1]
    fa_w_ew[:, :, 1:] = fa_e_ew[:, :, :-1]
    fa_w_ew[:, :, 0] = fa_e_ew[:, :, -1]

    return fa_e_we, fa_e_ew, fa_w_we, fa_w_ew


def to_boundaries_meridional(fa_n, p):
    """Define the horizontal fluxes over the north/south boundaries."""
    fa_n_boundary = np.zeros_like(fa_n)
    fa_n_boundary[:, 1:, :] = 0.5 * (fa_n[:, :-1, :] + fa_n[:, 1:, :])

    # find out where the positive and negative fluxes are
    fa_n_pos = np.ones_like(fa_n)
    fa_n_pos[fa_n_boundary < 0] = 0  # invalid value encountered in less
    fa_n_neg = fa_n_pos - 1

    # separate directions south-north (all positive numbers)
    fa_n_sn = fa_n_boundary * fa_n_pos
    fa_n_ns = fa_n_boundary * fa_n_neg

    # fluxes over the southern boundary
    fa_s_sn = np.nan * np.zeros_like(p)
    fa_s_ns = np.nan * np.zeros_like(p)
    fa_s_sn[:, :-1, :] = fa_n_sn[:, 1:, :]
    fa_s_ns[:, :-1, :] = fa_n_ns[:, 1:, :]

    return fa_n_sn, fa_n_ns, fa_s_sn, fa_s_ns


def backtrack(
    latitude,
    longitude,
    Kvf,
    region,
    Fa_east_upper,
    Fa_north_upper,
    Fa_east_lower,
    Fa_north_lower,
    Fa_Vert,
    E,
    P,
    s_upper,
    s_lower,
    s_track_upper_last,
    s_track_lower_last,
):
    P_region = region * P

    # Total moisture in the column
    W = s_upper + s_lower

    # separate the direction of the vertical flux and make it absolute
    Fa_upward = np.zeros(np.shape(Fa_Vert))
    Fa_upward[Fa_Vert <= 0] = Fa_Vert[Fa_Vert <= 0]
    Fa_downward = np.zeros(np.shape(Fa_Vert))
    Fa_downward[Fa_Vert >= 0] = Fa_Vert[Fa_Vert >= 0]
    Fa_upward = np.abs(Fa_upward)

    # include the vertical dispersion
    if Kvf != 0:
        Fa_upward = (1.0 + Kvf) * Fa_upward
        Fa_upward[Fa_Vert >= 0] = Fa_Vert[Fa_Vert >= 0] * Kvf
        Fa_downward = (1.0 + Kvf) * Fa_downward
        Fa_downward[Fa_Vert <= 0] = np.abs(Fa_Vert[Fa_Vert <= 0]) * Kvf

    # Define the horizontal fluxes over the boundaries
    Fa_east_upper_WE, Fa_east_upper_EW, Fa_west_upper_WE, Fa_west_upper_EW = to_boundaries_zonal(Fa_east_upper, P)
    Fa_east_lower_WE, Fa_east_lower_EW, Fa_west_lower_WE, Fa_west_lower_EW = to_boundaries_zonal(Fa_east_lower, P)
    Fa_north_upper_SN, Fa_north_upper_NS, Fa_south_upper_SN, Fa_south_upper_NS = to_boundaries_meridional(Fa_north_upper, P)
    Fa_north_lower_SN, Fa_north_lower_NS, Fa_south_lower_SN, Fa_south_lower_NS = to_boundaries_meridional(Fa_north_lower, P)

    # defining size of output
    s_track_lower = np.zeros_like(s_lower)
    s_track_upper = np.zeros_like(s_upper)

    # assign begin values of output == last (but first index) values of the previous time slot
    s_track_lower[-1, :, :] = s_track_lower_last
    s_track_upper[-1, :, :] = s_track_upper_last

    # defining sizes of tracked moisture
    s_track_after_Fa_lower = np.zeros_like(s_track_lower_last)
    s_track_after_Fa_P_east_lower = np.zeros_like(s_track_lower_last)
    s_track_east_lower = np.zeros_like(s_track_lower_last)
    s_track_west_lower = np.zeros_like(s_track_lower_last)
    s_track_north_lower = np.zeros_like(s_track_lower_last)
    s_track_south_lower = np.zeros_like(s_track_lower_last)
    s_track_after_Fa_upper = np.zeros_like(s_track_upper_last)
    s_track_after_Fa_P_east_upper = np.zeros_like(s_track_upper_last)
    s_track_east_upper = np.zeros_like(s_track_upper_last)
    s_track_west_upper = np.zeros_like(s_track_upper_last)
    s_track_north_upper = np.zeros_like(s_track_upper_last)
    s_track_south_upper = np.zeros_like(s_track_upper_last)

    # define sizes of total moisture
    s_east_lower = np.zeros_like(s_track_lower_last)
    s_west_lower = np.zeros_like(s_track_lower_last)
    s_north_lower = np.zeros_like(s_track_lower_last)
    s_south_lower = np.zeros_like(s_track_lower_last)
    s_east_upper = np.zeros_like(s_track_upper_last)
    s_west_upper = np.zeros_like(s_track_upper_last)
    s_north_upper = np.zeros_like(s_track_upper_last)
    s_south_upper = np.zeros_like(s_track_upper_last)

    # define variables that find out what happens to the water
    ntime = w_upper.shape[0]
    north_loss = np.zeros((ntime, 1, len(longitude)))
    south_loss = np.zeros((ntime, 1, len(longitude)))
    east_loss = np.zeros((ntime, 1, len(latitude)))
    west_loss = np.zeros((ntime, 1, len(latitude)))
    down_to_upper = np.zeros_like(P)
    top_to_lower = np.zeros_like(P)
    water_lost = np.zeros_like(P)

    # Sa calculation backward in time
    for t in reversed(range(ntime)):

        # Lower layer: define values of total moisture
        s_north_lower[0, 1:, :] = s_lower[t, :-1, :]
        s_south_lower[0, :-1, :] = s_lower[t, 1:, :]
        s_east_lower[0, :, :-1] = s_lower[t, :, 1:]
        s_west_lower[0, :, 1:] = s_lower[t, :, :-1]
        # TODO: fix div by zero issue for periodic_boundary
        # s_east_lower[0,:,-1] = s_lower[t,:,0]
        # s_west_lower[0,:,0] = s_lower[t,:,-1]

        # Lower layer: define values of tracked moisture of neighbouring grid cells
        s_track_north_lower[0, 1:, :] = s_track_lower[t, :-1, :]
        s_track_south_lower[0, :-1, :] = s_track_lower[t, 1:, :]
        s_track_east_lower[0, :, :-1] = s_track_lower[t, :, 1:]
        s_track_west_lower[0, :, 1:] = s_track_lower[t, :, :-1]
        if config["periodic_boundary"]:
            s_track_east_lower[0, :, -1] = s_track_lower[t, :, 0]
            s_track_west_lower[0, :, 0] = s_track_lower[t, :, -1]

        # Lower layer: calculate with moisture fluxes
        s_track_after_Fa_lower[0, 1:-1, 1:-1] = (
            s_track_lower[t, 1:-1, 1:-1]
            + Fa_east_lower_WE[t - 1, 1:-1, 1:-1] * (s_track_east_lower[0, 1:-1, 1:-1] / s_east_lower[0, 1:-1, 1:-1])
            - Fa_east_lower_EW[t - 1, 1:-1, 1:-1] * (s_track_lower[t, 1:-1, 1:-1] / s_lower[t, 1:-1, 1:-1])
            - Fa_west_lower_WE[t - 1, 1:-1, 1:-1] * (s_track_lower[t, 1:-1, 1:-1] / s_lower[t, 1:-1, 1:-1])
            + Fa_west_lower_EW[t - 1, 1:-1, 1:-1] * (s_track_west_lower[0, 1:-1, 1:-1] / s_west_lower[0, 1:-1, 1:-1])
            + Fa_north_lower_SN[t - 1, 1:-1, 1:-1] * (s_track_north_lower[0, 1:-1, 1:-1] / s_north_lower[0, 1:-1, 1:-1])
            - Fa_north_lower_NS[t - 1, 1:-1, 1:-1] * (s_track_lower[t, 1:-1, 1:-1] / s_lower[t, 1:-1, 1:-1])
            - Fa_south_lower_SN[t - 1, 1:-1, 1:-1] * (s_track_lower[t, 1:-1, 1:-1] / s_lower[t, 1:-1, 1:-1])
            + Fa_south_lower_NS[t - 1, 1:-1, 1:-1] * (s_track_south_lower[0, 1:-1, 1:-1] / s_south_lower[0, 1:-1, 1:-1])
            - Fa_downward[t - 1, 1:-1, 1:-1] * (s_track_lower[t, 1:-1, 1:-1] / s_lower[t, 1:-1, 1:-1])
            + Fa_upward[t - 1, 1:-1, 1:-1] * (s_track_upper[t, 1:-1, 1:-1] / s_upper[t, 1:-1, 1:-1])
        )

        # Upper layer: define values of total moisture
        s_north_upper[0, 1:, :] = s_upper[t, :-1, :]
        s_south_upper[0, :-1, :] = s_upper[t, 1:, :]
        s_east_upper[0, :, :-1] = s_upper[t, :, 1:]
        s_west_upper[0, :, 1:] = s_upper[t, :, :-1]
        # TODO: fix div by zero issue for periodic_boundary
        # s_east_upper[0,:,-1] = s_upper[t,:,0]
        # s_west_upper[0,:,0] = s_upper[t,:,-1]

        # Upper layer: define values of tracked moisture of neighbouring grid cells
        s_track_north_upper[0, 1:, :] = s_track_upper[t, :-1, :]
        s_track_south_upper[0, :-1, :] = s_track_upper[t, 1:, :]
        s_track_east_upper[0, :, :-1] = s_track_upper[t, :, 1:]
        s_track_west_upper[0, :, 1:] = s_track_upper[t, :, :-1]
        if config["periodic_boundary"]:
            s_track_east_upper[0, :, -1] = s_track_upper[t, :, 0]
            s_track_west_upper[0, :, 0] = s_track_upper[t, :, -1]

        # top: calculate with moisture fluxes
        s_track_after_Fa_upper[0, 1:-1, 1:-1] = (
            s_track_upper[t, 1:-1, 1:-1]
            + Fa_east_upper_WE[t - 1, 1:-1, 1:-1] * (s_track_east_upper[0, 1:-1, 1:-1] / s_east_upper[0, 1:-1, 1:-1])
            - Fa_east_upper_EW[t - 1, 1:-1, 1:-1] * (s_track_upper[t, 1:-1, 1:-1] / s_upper[t, 1:-1, 1:-1])
            - Fa_west_upper_WE[t - 1, 1:-1, 1:-1] * (s_track_upper[t, 1:-1, 1:-1] / s_upper[t, 1:-1, 1:-1])
            + Fa_west_upper_EW[t - 1, 1:-1, 1:-1] * (s_track_west_upper[0, 1:-1, 1:-1] / s_west_upper[0, 1:-1, 1:-1])
            + Fa_north_upper_SN[t - 1, 1:-1, 1:-1] * (s_track_north_upper[0, 1:-1, 1:-1] / s_north_upper[0, 1:-1, 1:-1])
            - Fa_north_upper_NS[t - 1, 1:-1, 1:-1] * (s_track_upper[t, 1:-1, 1:-1] / s_upper[t, 1:-1, 1:-1])
            - Fa_south_upper_SN[t - 1, 1:-1, 1:-1] * (s_track_upper[t, 1:-1, 1:-1] / s_upper[t, 1:-1, 1:-1])
            + Fa_south_upper_NS[t - 1, 1:-1, 1:-1] * (s_track_south_upper[0, 1:-1, 1:-1] / s_south_upper[0, 1:-1, 1:-1])
            + Fa_downward[t - 1, 1:-1, 1:-1] * (s_track_lower[t, 1:-1, 1:-1] / s_lower[t, 1:-1, 1:-1])
            - Fa_upward[t - 1, 1:-1, 1:-1] * (s_track_upper[t, 1:-1, 1:-1] / s_upper[t, 1:-1, 1:-1])
        )

        # losses to the north and south
        north_loss[t - 1, 0, :] = Fa_north_upper_NS[t - 1, 1, :] * (s_track_upper[t, 1, :] / s_upper[t, 1, :]) + Fa_north_lower_NS[t - 1, 1, :] * (s_track_lower[t, 1, :] / s_lower[t, 1, :])
        south_loss[t - 1, 0, :] = Fa_south_upper_SN[t - 1, -2, :] * (s_track_upper[t, -2, :] / s_upper[t, -2, :]) + Fa_south_lower_SN[t - 1, -2, :] * (s_track_lower[t, -2, :] / s_lower[t, -2, :])

        # losses to the east and west
        east_loss[t - 1, 0, :] = Fa_east_upper_EW[t - 1, :, -2] * (s_track_upper[t, :, -2] / s_upper[t, :, -2]) + Fa_east_lower_EW[t - 1, :, -2] * (s_track_lower[t, :, -2] / s_lower[t, :, -2])
        west_loss[t - 1, 0, :] = Fa_west_upper_WE[t - 1, :, 1] * (s_track_upper[t, :, 1] / s_upper[t, :, 1]) + Fa_west_lower_WE[t - 1, :, 1] * (s_track_lower[t, :, 1] / s_lower[t, :, 1])

        # down: add precipitation and subtract evaporation
        s_track_after_Fa_P_east_lower[0, 1:-1, 1:-1] = (
            s_track_after_Fa_lower[0, 1:-1, 1:-1]
            + P_region[t - 1, 1:-1, 1:-1] * (s_lower[t, 1:-1, 1:-1] / W[t, 1:-1, 1:-1])
            - E[t - 1, 1:-1, 1:-1] * (s_track_lower[t, 1:-1, 1:-1] / s_lower[t, 1:-1, 1:-1])
        )

        # top: add precipitation
        s_track_after_Fa_P_east_upper[0, 1:-1, 1:-1] = s_track_after_Fa_upper[0, 1:-1, 1:-1] + P_region[t - 1, 1:-1, 1:-1] * (s_upper[t, 1:-1, 1:-1] / W[t, 1:-1, 1:-1])

        # down and top: redistribute unaccounted water that is otherwise lost from the sytem
        down_to_upper[t - 1, :, :] = np.reshape(np.maximum(0, np.reshape(s_track_after_Fa_P_east_lower, (np.size(s_track_after_Fa_P_east_lower))) - np.reshape(s_lower[t - 1, :, :], (np.size(s_lower[t - 1, :, :]))),),(len(latitude), len(longitude)),)  # Imme: invalid value encountered in maximum
        top_to_lower[t - 1, :, :] = np.reshape(np.maximum(0, np.reshape(s_track_after_Fa_P_east_upper, (np.size(s_track_after_Fa_P_east_upper))) - np.reshape(s_upper[t - 1, :, :], (np.size(s_upper[t - 1, :, :]))),),(len(latitude), len(longitude)),)  # Imme: invalid value encountered in maximum
        s_track_lower[t - 1, 1:-1, 1:-1] = (s_track_after_Fa_P_east_lower[0, 1:-1, 1:-1] - down_to_upper[t - 1, 1:-1, 1:-1] + top_to_lower[t - 1, 1:-1, 1:-1])
        s_track_upper[t - 1, 1:-1, 1:-1] = (s_track_after_Fa_P_east_upper[0, 1:-1, 1:-1] - top_to_lower[t - 1, 1:-1, 1:-1] + down_to_upper[t - 1, 1:-1, 1:-1])

    return (
        s_track_upper,
        s_track_lower,
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
