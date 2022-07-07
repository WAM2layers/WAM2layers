import yaml
import xarray as xr

import numpy as np
from pathlib import Path
import pandas as pd

from preprocessing import get_grid_info, stabilize_fluxes

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


def resample(ds, target_freq):
    """Increase time resolution; states at midpoints, fluxes at the edges."""
    target_seconds = pd.Timedelta(target_freq).total_seconds()
    current_seconds = ds.time.diff('time').dt.seconds.values[0]
    resample_ratio = current_seconds / target_seconds

    time = ds.time.values
    newtime_states = pd.date_range(time[0], time[-1], freq=target_freq)
    newtime_fluxes = newtime_states[:-1] + pd.Timedelta(target_freq) / 2

    states = ds[['s_upper', 's_lower']].interp(time=newtime_states)
    fluxes = ds[['fx_upper', 'fx_lower', 'fy_upper', 'fy_lower']].interp(time=newtime_fluxes)
    surface = ds[['precip', 'evap']].reindex(time=newtime_fluxes, method="bfill") / resample_ratio
    return fluxes.merge(surface), states


def change_units(fluxes, target_freq):
    """Change units to m3.

    Multiply by edge length to get flux in m3
    Multiply by time to get accumulation instead of flux
    Divide by density of water to go from kg to m3
    """
    density = 1000  # [kg/m3]
    a, ly, lx = get_grid_info(fluxes)

    total_seconds = pd.Timedelta(target_freq).total_seconds()
    fluxes.fx_upper *= total_seconds / density * ly
    fluxes.fx_lower *= total_seconds / density * ly
    fluxes.fy_upper *= total_seconds / density * lx[None, :, None]
    fluxes.fy_lower *= total_seconds / density * lx[None, :, None]
    fluxes.evap *= a
    fluxes.precip *= a

    for variable in fluxes.data_vars:
        fluxes[variable].assign_units("m**3")


def to_edges_zonal(fx):
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


def to_edges_meridional(fy):
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


def split_vertical_flux(Kvf, fv):
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

    return f_downward, f_upward


def backtrack(
    preprocessed_data,
    s_track_upper,
    s_track_lower,
    region,
    kvf,
):

    # Unpack preprocessed data  # TODO: change names in preproc
    fx_upper = preprocessed_data["fx_upper"].values
    fy_upper = preprocessed_data["fy_upper"].values
    fx_lower = preprocessed_data["fx_lower"].values
    fy_lower = preprocessed_data["fy_lower"].values
    evap = preprocessed_data["evap"].values
    precip = preprocessed_data["precip"].values
    s_upper = preprocessed_data["s_upper"].values
    s_lower = preprocessed_data["s_lower"].values
    f_vert = preprocessed_data["f_vert"].values

    # Allocate arrays for daily accumulations
    ntime, nlat, nlon = s_upper.shape

    s_track_upper_agg = np.zeros((nlat, nlon))
    s_track_lower_agg = np.zeros((nlat, nlon))
    e_track = np.zeros((nlat, nlon))

    north_loss = south_loss = np.zeros(nlon)
    east_loss = west_loss = np.zeros(nlat)

    # Sa calculation backward in time
    ntime = len(s_lower)
    for t in reversed(range(ntime)):
        P_region = region * precip[t-1]
        s_total = s_upper[t] + s_lower[t]

        # separate the direction of the vertical flux and make it absolute
        f_downward, f_upward = split_vertical_flux(kvf, f_vert[t-1])

        # Determine horizontal fluxes over the grid-cell boundaries
        fx_e_lower_we, fx_e_lower_ew, fx_w_lower_we, fx_w_lower_ew = to_edges_zonal(fx_lower[t-1])
        fx_e_upper_we, fx_e_upper_ew, fx_w_upper_we, fx_w_upper_ew = to_edges_zonal(fx_upper[t-1])
        fy_n_lower_sn, fy_n_lower_ns, fy_s_lower_sn, fy_s_lower_ns = to_edges_meridional(fy_lower[t-1])
        fy_n_upper_sn, fy_n_upper_ns, fy_s_upper_sn, fy_s_upper_ns = to_edges_meridional(fy_upper[t-1])

        # Short name for often used expressions
        s_track_relative_lower = s_track_lower / s_lower[t]  # fraction of tracked relative to total moisture
        s_track_relative_upper = s_track_upper / s_upper[t]
        inner = np.s_[1:-1, 1:-1]

        # Actual tracking (note: backtracking, all terms have been negated)
        s_track_lower[inner] += (
            + fx_e_lower_we * shift_east(s_track_relative_lower)
            + fx_w_lower_ew * shift_west(s_track_relative_lower)
            + fy_n_lower_sn * shift_north(s_track_relative_lower)
            + fy_s_lower_ns * shift_south(s_track_relative_lower)
            + f_upward * s_track_relative_upper
            - f_downward * s_track_relative_lower
            - fy_s_lower_sn * s_track_relative_lower
            - fy_n_lower_ns * s_track_relative_lower
            - fx_e_lower_ew * s_track_relative_lower
            - fx_w_lower_we * s_track_relative_lower
            + P_region * (s_lower[t] / s_total)
            - evap[t-1] * s_track_relative_lower
        )[inner]

        s_track_upper[inner] += (
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
        lower_to_upper = np.maximum(0, s_track_lower) - s_lower[t - 1]
        upper_to_lower = np.maximum(0, s_track_upper) - s_upper[t - 1]
        s_track_lower[inner] = (s_track_lower - lower_to_upper + upper_to_lower)[inner]
        s_track_upper[inner] = (s_track_upper - upper_to_lower + lower_to_upper)[inner]

        # compute tracked evaporation
        e_track += evap[t-1] * (s_track_lower / s_lower[t])

        # losses to the north and south
        north_loss += (fy_n_upper_ns * s_track_relative_upper
                       + fy_n_lower_ns * s_track_relative_lower)[1, :]

        south_loss += (fy_s_upper_sn * s_track_relative_upper
                       + fy_s_lower_sn * s_track_relative_lower)[-2, :]

        east_loss += (fx_e_upper_ew * s_track_relative_upper
                      + fx_e_lower_ew * s_track_relative_lower)[:, -2]

        west_loss += (fx_w_upper_we * s_track_relative_upper
                      + fx_w_lower_we * s_track_relative_lower)[:, 1]

        # Aggregate daily accumulations
        s_track_lower_agg += s_track_lower
        s_track_upper_agg += s_track_upper

    # Pack processed data into new dataset
    ds = xr.Dataset(
        {
            "s_track_upper_last": (["lat", "lon"], s_track_upper),  # Keep last state for a restart
            "s_track_lower_last": (["lat", "lon"], s_track_lower),
            "s_track_upper": (["lat", "lon"], s_track_upper_agg),
            "s_track_lower": (["lat", "lon"], s_track_lower_agg),
            "e_track": (["lat", "lon"], e_track),
            "north_loss": (["lon"], north_loss),
            "south_loss": (["lon"], south_loss),
            "east_loss": (["lat"], east_loss),
            "west_loss": (["lat",], west_loss),
        }
    )
    return (
        s_track_upper,
        s_track_lower,
        ds
    )


region = xr.open_dataset(config['region']).isel(time=0).lsm.rename(latitude='lat', longitude='lon').values

for i, date in enumerate(reversed(datelist)):
    print(date)

    if i == 0:
        if config["restart"]:
            # Reload last state from existing output
            ds = xr.open_dataset(output_path(date + pd.Timedelta(days=1)))
            s_track_upper = ds.s_track_upper.values
            s_track_lower = ds.s_track_lower.values
        else:
            # Allocate empty arrays based on shape of input data
            ds = xr.open_dataset(input_path(datelist[0]))
            s_track_upper = np.zeros_like(ds.s_upper[0])
            s_track_lower = np.zeros_like(ds.s_upper[0])

    preprocessed_data = xr.open_dataset(input_path(date))

    # Resample to (higher) target frequency
    # After this, the fluxes will be "in between" the states
    fluxes, states = resample(preprocessed_data, config['target_frequency'])

    # Convert flux data to volumes
    change_units(fluxes, config["target_frequency"])

    # TODO: Try stacking upper and lower into a single array with new dim
    # Stabilize horizontal fluxes # TODO reinsert in datasets
    fx_upper, fy_upper = stabilize_fluxes(fluxes.fx_upper, fluxes.fy_upper, states.s_upper)
    fx_lower, fy_lower = stabilize_fluxes(fx_lower, fy_lower, s_lower)

    # Determine the vertical moisture flux
    f_vert = get_vertical_transport(ds_fluxes, ds_accumulations, ds_states,
                                    config["periodic_boundary"], config["kvf"])

    (s_track_upper, s_track_lower, processed_data) = backtrack(
        preprocessed_data,
        s_track_upper,
        s_track_lower,
        region,
        config['kvf'],
    )

    # Write output to file
    # TODO: add (and cleanup) coordinates and units
    processed_data.to_netcdf(output_path(date))
