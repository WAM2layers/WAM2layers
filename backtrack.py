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


def to_boundaries_zonal(fx):
    """Define the horizontal fluxes over the east/west boundaries."""
    fxh = np.zeros_like(fx)
    fxh[:, :-1] = 0.5 * (fx[:, :-1] + fx[:, 1:])
    if config["periodic_boundary"]:
        fxh[:, -1] = 0.5 * (fx[:, -1] + fx[:, 0])

    # find out where the positive and negative fluxes are
    fx_pos = np.ones_like(fx)
    fx_pos[fxh < 0] = 0
    fx_neg = fx_pos - 1

    # separate directions west-east (all positive numbers)
    fx_e_we = fxh * fx_pos  # eastern edge, outgoing to west
    fx_e_ew = fxh * fx_neg  # eastern edge, incoming from west

    # fluxes over the western boundary
    fx_w_we = shift_west(fx_e_we)
    fx_w_ew = shift_west(fx_e_ew)

    return fx_e_we, fx_e_ew, fx_w_we, fx_w_ew


def to_boundaries_meridional(fy):
    """Define the horizontal fluxes over the north/south boundaries."""
    fy_boundary = np.zeros_like(fy)
    fy_boundary[1:, :] = 0.5 * (fy[:-1, :] + fy[1:, :])

    # find out where the positive and negative fluxes are
    fy_pos = np.ones_like(fy)
    fy_pos[fy_boundary < 0] = 0  # invalid value encountered in less
    fy_neg = fy_pos - 1

    # separate directions south-north (all positive numbers)
    fy_n_sn = fy_boundary * fy_pos
    fy_n_ns = fy_boundary * fy_neg

    # fluxes over the southern boundary
    # TODO: shouldn't this be shift_south?
    fy_s_sn = shift_north(fy_n_sn)
    fy_s_ns = shift_north(fy_n_ns)

    return fy_n_sn, fy_n_ns, fy_s_sn, fy_s_ns


def shift_north(array):
    # Note: edges are reinserted at other end; but they're not used anyway
    return np.roll(array, -1, axis=-2)


def shift_south(array):
    # Note: edges are reinserted at other end; but they're not used anyway
    return np.roll(array, 1, axis=-2)


def shift_east(array):
    # Note: edges are reinserted at other end; but they're not used anyway
    return np.roll(array, -1, axis=-1)


def shift_west(array):
    # Note: edges are reinserted at other end; but they're not used anyway
    return np.roll(array, 1, axis=-1)


def backtrack(
    Kvf,
    region,
    fx_upper,
    fy_upper,
    fx_lower,
    fy_lower,
    f_vert,
    E,
    P,
    s_upper,
    s_lower,
    s_track_upper_last,
    s_track_lower_last,
):

    # assign begin values of output == last (but first index) values of the previous time slot
    s_track_lower = s_track_lower_last.squeeze()
    s_track_upper = s_track_upper_last.squeeze()

    # Sa calculation backward in time
    ntime = len(s_lower)
    for t in reversed(range(ntime)):
        P_region = region * P[t-1]
        s_total = s_upper[t] + s_lower[t]

        # separate the direction of the vertical flux and make it absolute
        fv = f_vert[t-1]
        f_downward = f_upward = np.zeros_like(fv)
        f_downward[fv >= 0] = fv[fv >= 0]
        f_upward[fv <= 0] = fv[fv <= 0]
        f_upward = np.abs(f_upward)

        # include the vertical dispersion
        if Kvf != 0:
            f_upward = (1.0 + Kvf) * f_upward
            f_upward[fv >= 0] = fv[fv >= 0] * Kvf
            f_downward = (1.0 + Kvf) * f_downward
            f_downward[fv <= 0] = np.abs(fv[fv <= 0]) * Kvf

        # Define the horizontal fluxes over the boundaries
        fx_e_upper_we, fx_e_upper_ew, fx_w_upper_we, fx_w_upper_ew = to_boundaries_zonal(fx_upper[t-1])
        fx_e_lower_we, fx_e_lower_ew, fx_w_lower_we, fx_w_lower_ew = to_boundaries_zonal(fx_lower[t-1])
        fy_n_upper_sn, fy_n_upper_ns, fy_s_upper_sn, fy_s_upper_ns = to_boundaries_meridional(fy_upper[t-1])
        fy_n_lower_sn, fy_n_lower_ns, fy_s_lower_sn, fy_s_lower_ns = to_boundaries_meridional(fy_lower[t-1])

        # Short name for often used expressions
        s_track_relative_lower = s_track_lower / s_lower[t]  # fraction of tracked relative to total moisture
        s_track_relative_upper = s_track_upper / s_upper[t]
        inner = np.s_[1:-1, 1:-1]

        # Actual tracking
        s_track_lower[inner] = (
            s_track_lower  # current state
            # TODO: Looks like I'm messing up my interpretation of incoming/outgoing here...
            # I think fx_e_we should be an outgoing flux, but I'm not sure
            + fx_e_lower_we * shift_east(s_track_relative_lower)  # in from the east
            + fx_w_lower_ew * shift_west(s_track_relative_lower)  # in from the west
            + fy_n_lower_sn * shift_north(s_track_relative_lower)  # in from the north
            + fy_s_lower_ns * shift_south(s_track_relative_lower)  # in from the south
            + f_upward * s_track_relative_upper  # in from upper layer
            - f_downward * s_track_relative_lower  # out to the upper layer
            - fy_s_lower_sn * s_track_relative_lower  # out to the south
            - fy_n_lower_ns * s_track_relative_lower  # out to the north
            - fx_e_lower_ew * s_track_relative_lower  # out to the west
            - fx_w_lower_we * s_track_relative_lower  # out to the east
            + P_region * (s_lower[t] / s_total)  # gained from precipitation?
            - E[t-1] * s_track_relative_lower  # lost to evaporation?
        )[inner]

        s_track_upper[inner] = (
            s_track_upper
            + fx_e_upper_we * shift_east(s_track_relative_upper)
            + fx_w_upper_ew * shift_west(s_track_relative_upper)
            + fy_n_upper_sn * shift_north(s_track_relative_upper)
            + fy_s_upper_ns * shift_south(s_track_relative_upper)
            + f_downward * s_track_relative_lower
            - f_upward * s_track_relative_upper
            - fy_s_upper_sn * s_track_relative_upper
            - fy_n_upper_ns * s_track_relative_upper
            - fx_w_upper_we * s_track_relative_upper
            - fx_e_upper_ew * s_track_relative_upper
            + P_region * (s_upper[t] / s_total)
        )[inner]

        # down and top: redistribute unaccounted water that is otherwise lost from the sytem
        down_to_upper = np.maximum(0, s_track_lower) - s_lower[t - 1, :, :]
        top_to_lower = np.maximum(0, s_track_upper) - s_upper[t - 1, :, :]
        s_track_lower[inner] = (s_track_lower - down_to_upper + top_to_lower)[inner]
        s_track_upper[inner] = (s_track_upper - top_to_lower + down_to_upper)[inner]

        # losses to the north and south
        north_loss = (fy_n_upper_ns * s_track_relative_upper
                      + fy_n_lower_ns * s_track_relative_lower)[1, :]

        south_loss = (fy_s_upper_sn * s_track_relative_upper
                      + fy_s_lower_sn * s_track_relative_lower)[-2, :]

        east_loss = (fx_e_upper_ew * s_track_relative_upper
                     + fx_e_lower_ew * s_track_relative_lower)[:, -2]

        west_loss = (fx_w_upper_we * s_track_relative_upper
                     + fx_w_lower_we * s_track_relative_lower)[:, 1]

    return (
        s_track_upper,
        s_track_lower,
        north_loss,
        south_loss,
        east_loss,
        west_loss,
        down_to_upper,
        top_to_lower,
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
        down_to_upper, top_to_lower) = backtrack(
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
        }
    ).to_netcdf(output_path(date))
