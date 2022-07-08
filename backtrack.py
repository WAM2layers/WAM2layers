import yaml
import xarray as xr

import numpy as np
from pathlib import Path
import pandas as pd

from preprocessing import get_grid_info, get_vertical_transport, stabilize_fluxes

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

    states = ds.s.interp(time=newtime_states)
    fluxes = ds[['fx', 'fy']].interp(time=newtime_fluxes)
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
    fluxes['fx'] *= total_seconds / density * ly
    fluxes['fy'] *= total_seconds / density * lx[None, :, None]
    fluxes['evap'] *= a[None, :, None]
    fluxes['precip'] *= a[None, :, None]

    for variable in fluxes.data_vars:
        fluxes[variable].assign_attrs(units="m**3")


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
    s_track,
    preprocessed_data,
    region,
    kvf,
):

    # Unpack preprocessed data  # TODO: change names in preproc
    s = preprocessed_data["s"].values
    fx = preprocessed_data["fx"].values
    fy = preprocessed_data["fy"].values
    fz = preprocessed_data["fz"].values
    evap = preprocessed_data["evap"].values
    precip = preprocessed_data["precip"].values

    # Allocate arrays for daily accumulations
    _, ntime, nlat, nlon = s.shape

    s_track_mean = np.zeros((2, nlat, nlon))
    e_track = np.zeros((nlat, nlon))

    north_loss = south_loss = np.zeros(nlon)
    east_loss = west_loss = np.zeros(nlat)

    # Sa calculation backward in time
    for t in reversed(range(ntime)):
        P_region = region * precip[t-1]
        s_total = s.sum(axis=0)[t]

        # separate the direction of the vertical flux and make it absolute
        f_downward, f_upward = split_vertical_flux(kvf, fz[t-1])

        # Determine horizontal fluxes over the grid-cell boundaries
        fx_e_we, fx_e_ew, fx_w_we, fx_w_ew = to_edges_zonal(fx[:, t-1])
        fy_n_sn, fy_n_ns, fy_s_sn, fy_s_ns = to_edges_meridional(fy[:, t-1])

        # Short name for often used expressions
        s_track_relative = s_track / s[:, t, ...]  # fraction of tracked relative to total moisture
        inner = np.s_[:, 1:-1, 1:-1]

        # Actual tracking (note: backtracking, all terms have been negated)
        s_track[inner] += (
            + fx_e_we * shift_east(s_track_relative)
            + fx_w_ew * shift_west(s_track_relative)
            + fy_n_sn * shift_north(s_track_relative)
            + fy_s_ns * shift_south(s_track_relative)
            - fy_s_sn * s_track_relative
            - fy_n_ns * s_track_relative
            - fx_e_ew * s_track_relative
            - fx_w_we * s_track_relative
            + P_region * (s[:, t, ...] / s_total)
        )[inner]

        # Specific to lower layer
        s_track[0, 1:-1, 1:-1] += (
            + f_upward * s_track_relative[1]
            - f_downward * s_track_relative[0]
            - evap[t-1] * s_track_relative[0]
        )[1:-1, 1:-1]

        # Specific to upper layer
        s_track[1, 1:-1, 1:-1] += (
            - f_upward * s_track_relative[1]
            + f_downward * s_track_relative[0]
        )[1:-1, 1:-1]

        # down and top: redistribute unaccounted water that is otherwise lost from the sytem
        # TODO: friday afternoon, might have messed this up.
        debet = np.maximum(0, s_track - s[:, t - 1])
        s_track[inner] = (s_track - debet[::-1])[inner]

        # compute tracked evaporation
        e_track += evap[t-1] * (s_track[0] / s[0][t])

        # losses to the north and south
        north_loss += (fy_n_ns * s_track_relative).sum(axis=0)[1, :]
        south_loss += (fy_s_sn * s_track_relative).sum(axis=0)[-2, :]
        east_loss += (fx_e_ew * s_track_relative).sum(axis=0)[:, -2]
        west_loss += (fx_w_we * s_track_relative).sum(axis=0)[:, 1]

        # Aggregate to eventually get a daily mean accumulations
        s_track_mean += s_track / ntime

    # Pack processed data into new dataset
    ds = xr.Dataset(
        {
            "s_track_restart": (["lev", "lat", "lon"], s_track),  # Keep last state for a restart
            "s_track": (["lev", "lat", "lon"], s_track_mean),
            "e_track": (["lat", "lon"], e_track),
            "north_loss": (["lon"], north_loss),
            "south_loss": (["lon"], south_loss),
            "east_loss": (["lat"], east_loss),
            "west_loss": (["lat",], west_loss),
        }
    )
    return s_track, ds


region = xr.open_dataset(config['region']).isel(time=0).lsm.rename(latitude='lat', longitude='lon').values

for i, date in enumerate(reversed(datelist)):
    print(date)

    if i == 0:
        if config["restart"]:
            # Reload last state from existing output
            ds = xr.open_dataset(output_path(date + pd.Timedelta(days=1)))
            s_track = ds.s_track_restart.values
        else:
            # Allocate empty arrays based on shape of input data
            ds = xr.open_dataset(input_path(datelist[0]))
            s_track = np.zeros_like(ds.s[:, 0, ...])

    preprocessed_data = xr.open_dataset(input_path(date))

    # Resample to (higher) target frequency
    # After this, the fluxes will be "in between" the states
    fluxes, s = resample(preprocessed_data, config['target_frequency'])

    # Convert flux data to volumes
    change_units(fluxes, config["target_frequency"])

    # Stabilize horizontal fluxes # TODO reinsert in datasets
    fx, fy = stabilize_fluxes(fluxes.fx, fluxes.fy, s)
    evap = fluxes.evap
    precip = fluxes.precip

    # Determine the vertical moisture flux
    kvf = config['kvf']
    pb = config['periodic_boundary']

    # Calculate fluxes and add to preprocessed data
    preprocessed_data['fz'] = get_vertical_transport(s, fx, fy, evap, precip, kvf, pb)

    s_track, processed_data = backtrack(s_track, preprocessed_data, region, config['kvf'])

    # Write output to file
    # TODO: add (and cleanup) coordinates and units
    processed_data.to_netcdf(output_path(date))
