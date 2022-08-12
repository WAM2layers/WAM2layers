from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import yaml

from visualization import make_diagnostic_figures
from preprocessing import get_grid_info

# Read case configuration
with open("cases/era5_2021.yaml") as f:
    config = yaml.safe_load(f)


datelist = pd.date_range(
    start=config["track_start_date"],
    end=config["track_end_date"],
    freq="d",
    inclusive="left",
)

input_dir = Path(config["preprocessed_data_folder"]).expanduser()
output_dir = Path(config["output_folder"]).expanduser() / "backtrack"

# Check if input dir exists
if not input_dir.exists():
    raise ValueError(
        "Please create the preprocessed_data_folder before running the script"
    )

# Create output dir if it doesn't exist yet
if not output_dir.exists():
    output_dir.mkdir(parents=True)


def time_in_range(start, end, current):
    """Returns whether current is in the range [start, end]"""
    return start <= current <= end


def input_path(date):
    return f"{input_dir}/{date.strftime('%Y-%m-%d')}_fluxes_storages.nc"


def output_path(date):
    return f"{output_dir}/{date.strftime('%Y-%m-%d')}_s_track.nc"


def to_edges_zonal(fx, periodic_boundary=False):
    """Define the horizontal fluxes over the east/west boundaries."""
    fe = np.zeros_like(fx)
    fe[:, :-1] = 0.5 * (fx[:, :-1] + fx[:, 1:])
    if periodic_boundary:
        fe[:, -1] = 0.5 * (fx[:, -1] + fx[:, 0])

    # find out where the positive and negative fluxes are
    f_pos = np.ones_like(fx)
    f_pos[fe < 0] = 0
    f_neg = f_pos - 1

    # separate directions west-east (all positive numbers)
    fe_we = fe * f_pos
    fe_ew = fe * f_neg

    # fluxes over the western boundary
    fw_we = look_west(fe_we)
    fw_ew = look_west(fe_ew)

    return fe_we, fe_ew, fw_we, fw_ew


def to_edges_meridional(fy):
    """Define the horizontal fluxes over the north/south boundaries."""
    fn = np.zeros_like(fy)
    fn[1:, :] = 0.5 * (fy[:-1, :] + fy[1:, :])

    # find out where the positive and negative fluxes are
    fn_pos = np.ones_like(fn)
    fn_pos[fn < 0] = 0  # invalid value encountered in less
    fn_neg = fn_pos - 1

    # separate directions south-north (all positive numbers)
    fn_sn = fn * fn_pos
    fn_ns = fn * fn_neg

    # fluxes over the southern boundary
    fs_sn = look_south(fn_sn)
    fs_ns = look_south(fn_ns)

    return fn_sn, fn_ns, fs_sn, fs_ns


def look_north(array):
    # Note: edges are reinserted at other end; but they're not used anyway
    return np.roll(array, 1, axis=-2)


def look_south(array):
    # Note: edges are reinserted at other end; but they're not used anyway
    return np.roll(array, -1, axis=-2)


def look_east(array):
    # Note: edges are reinserted at other end; but they're not used anyway
    return np.roll(array, -1, axis=-1)


def look_west(array):
    # Note: edges are reinserted at other end; but they're not used anyway
    return np.roll(array, 1, axis=-1)


def split_vertical_flux(Kvf, fv):
    f_downward = np.zeros_like(fv)
    f_upward = np.zeros_like(fv)
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
    fluxes["fx_upper"] *= total_seconds / density * ly
    fluxes["fx_lower"] *= total_seconds / density * ly
    fluxes["fy_upper"] *= total_seconds / density * lx[None, :, None]
    fluxes["fy_lower"] *= total_seconds / density * lx[None, :, None]
    fluxes["evap"] *= a[None, :, None]
    fluxes["precip"] *= a[None, :, None]

    for variable in fluxes.data_vars:
        fluxes[variable] = fluxes[variable].assign_attrs(units="m**3")


def stabilize_fluxes(fluxes, states):
    """Stabilize the outfluxes / influxes.

    During the reduced timestep the water cannot move further than 1/x * the
    gridcell, In other words at least x * the reduced timestep is needed to
    cross a gridcell.
    """
    for level in ["upper", "lower"]:
        fx = fluxes["fx_" + level]
        fy = fluxes["fy_" + level]
        s = states["s_" + level]

        fx_abs = np.abs(fx)
        fy_abs = np.abs(fy)
        ft_abs = fx_abs + fy_abs

        fx_corrected = 1/2 * fx_abs / ft_abs * s[:-1, :, :].values
        fx_stable = np.minimum(fx_abs, fx_corrected)

        fy_corrected = 1/2 * fy_abs / ft_abs * s[:-1, :, :].values
        fy_stable = np.minimum(fy_abs, fy_corrected)

        # Get rid of any nan values
        fx_stable.fillna(0)
        fy_stable.fillna(0)

        # Re-instate the sign
        fluxes["fx_"+ level] = np.sign(fx) * fx_stable
        fluxes["fy_"+ level] = np.sign(fy) * fy_stable


def convergence(fx, fy):
    # Note: latitude decreasing, hence positive fy gradient is convergence
    return np.gradient(fy, axis=-2) - np.gradient(fx, axis=-1)


def calculate_fv(fluxes, states, kvf, periodic):
    """Calculate the vertical fluxes.

    Note: fluxes are given at temporal midpoints between states.
    """
    s_total = states.s_upper + states.s_lower
    s_rel_upper = (states.s_upper / s_total).interp(time=fluxes.time)
    s_rel_lower = (states.s_lower / s_total).interp(time=fluxes.time)

    tendency_upper = convergence(fluxes.fx_upper, fluxes.fy_upper) - fluxes.precip.values * s_rel_upper
    tendency_lower = convergence(fluxes.fx_upper, fluxes.fy_upper) - fluxes.precip.values * s_rel_lower + fluxes.evap

    residual_upper = states.s_upper.diff("time").values - tendency_upper
    residual_lower = states.s_lower.diff("time").values - tendency_lower

    # compute the resulting vertical moisture flux; the vertical velocity so
    # that the new residual_lower/s_lower = residual_upper/s_upper (positive downward)
    fv = s_rel_lower * (residual_upper + residual_lower) - residual_lower

    # stabilize the outfluxes / influxes; during the reduced timestep the
    # vertical flux can maximally empty/fill 1/x of the top or down storage
    stab = 1.0 / (kvf + 1.0)
    flux_limit = np.minimum(states.s_upper, states.s_lower).interp(time=fluxes.time)
    fv_stable = np.minimum(np.abs(fv), stab * flux_limit)

    # Reinstate the sign
    return np.sign(fv) * fv_stable


def backtrack(
    date,
    fluxes,
    states,
    s_track_upper,
    s_track_lower,
    region,
    kvf,
):

    # Unpack preprocessed data
    fx_upper = fluxes["fx_upper"].values
    fy_upper = fluxes["fy_upper"].values
    fx_lower = fluxes["fx_lower"].values
    fy_lower = fluxes["fy_lower"].values
    evap = fluxes["evap"].values
    precip = fluxes["precip"].values
    f_vert = fluxes["f_vert"].values
    s_upper = states["s_upper"].values
    s_lower = states["s_lower"].values

    # Allocate arrays for daily accumulations
    ntime, nlat, nlon = fx_upper.shape

    s_track_upper_mean = np.zeros((nlat, nlon))
    s_track_lower_mean = np.zeros((nlat, nlon))
    e_track = np.zeros((nlat, nlon))

    north_loss = np.zeros(nlon)
    south_loss = np.zeros(nlon)
    east_loss = np.zeros(nlat)
    west_loss = np.zeros(nlat)

    # Only track the precipitation at certain dates
    if (
        time_in_range(
            config["event_start_date"],
            config["event_end_date"],
            date.strftime("%Y%m%d"),
        )
        == False
    ):
        precip = precip * 0

    # Sa calculation backward in time
    for t in reversed(range(ntime)):
        # die allerlaatste stap -1 wil je toch niet?????
        P_region = region * precip[t]  # WHY t-1 ?
        s_total = s_upper[t+1] + s_lower[t+1]

        # separate the direction of the vertical flux and make it absolute
        f_downward, f_upward = split_vertical_flux(kvf, f_vert[t])

        # Determine horizontal fluxes over the grid-cell boundaries
        f_e_lower_we, f_e_lower_ew, f_w_lower_we, f_w_lower_ew = to_edges_zonal(
            fx_lower[t]
        )
        f_e_upper_we, f_e_upper_ew, f_w_upper_we, f_w_upper_ew = to_edges_zonal(
            fx_upper[t]
        )

        (
            fy_n_lower_sn,
            fy_n_lower_ns,
            fy_s_lower_sn,
            fy_s_lower_ns,
        ) = to_edges_meridional(fy_lower[t])
        (
            fy_n_upper_sn,
            fy_n_upper_ns,
            fy_s_upper_sn,
            fy_s_upper_ns,
        ) = to_edges_meridional(fy_upper[t])

        # Short name for often used expressions
        s_track_relative_lower = (
            s_track_lower / s_lower[t+1]
        )  # fraction of tracked relative to total moisture
        s_track_relative_upper = s_track_upper / s_upper[t+1]
        inner = np.s_[1:-1, 1:-1]

        # Actual tracking (note: backtracking, all terms have been negated)
        s_track_lower[inner] += (
            + f_e_lower_we * look_east(s_track_relative_lower)
            + f_w_lower_ew * look_west(s_track_relative_lower)
            + fy_n_lower_sn * look_north(s_track_relative_lower)
            + fy_s_lower_ns * look_south(s_track_relative_lower)
            + f_upward * s_track_relative_upper
            - f_downward * s_track_relative_lower
            - fy_s_lower_sn * s_track_relative_lower
            - fy_n_lower_ns * s_track_relative_lower
            - f_e_lower_ew * s_track_relative_lower
            - f_w_lower_we * s_track_relative_lower
            + P_region * (s_lower[t+1] / s_total)
            - evap[t] * s_track_relative_lower
        )[inner]

        s_track_upper[inner] += (
            + f_e_upper_we * look_east(s_track_relative_upper)
            + f_w_upper_ew * look_west(s_track_relative_upper)
            + fy_n_upper_sn * look_north(s_track_relative_upper)
            + fy_s_upper_ns * look_south(s_track_relative_upper)
            + f_downward * s_track_relative_lower
            - f_upward * s_track_relative_upper
            - fy_s_upper_sn * s_track_relative_upper
            - fy_n_upper_ns * s_track_relative_upper
            - f_w_upper_we * s_track_relative_upper
            - f_e_upper_ew * s_track_relative_upper
            + P_region * (s_upper[t+1] / s_total)
        )[inner]

        # down and top: redistribute unaccounted water that is otherwise lost from the sytem
        lower_to_upper = np.maximum(0, s_track_lower - s_lower[t])
        upper_to_lower = np.maximum(0, s_track_upper - s_upper[t])
        s_track_lower[inner] = (s_track_lower - lower_to_upper + upper_to_lower)[inner]
        s_track_upper[inner] = (s_track_upper - upper_to_lower + lower_to_upper)[inner]

        # compute tracked evaporation
        e_track += evap[t] * (s_track_lower / s_lower[t+1])

        # losses to the north and south
        north_loss += (
            fy_n_upper_ns * s_track_relative_upper
            + fy_n_lower_ns * s_track_relative_lower
        )[1, :]

        south_loss += (
            fy_s_upper_sn * s_track_relative_upper
            + fy_s_lower_sn * s_track_relative_lower
        )[-2, :]

        east_loss += (
            f_e_upper_ew * s_track_relative_upper
            + f_e_lower_ew * s_track_relative_lower
        )[:, -2]

        west_loss += (
            f_w_upper_we * s_track_relative_upper
            + f_w_lower_we * s_track_relative_lower
        )[:, 1]

        # Aggregate daily accumulations for calculating the daily means
        s_track_lower_mean += s_track_lower / ntime
        s_track_upper_mean += s_track_upper / ntime

    make_diagnostic_figures(
        date,
        region,
        fx_upper,
        fy_upper,
        fx_lower,
        fy_lower,
        precip,
        s_track_upper_mean,
        s_track_lower_mean,
        e_track,
    )

    # Pack processed data into new dataset
    ds = xr.Dataset(
        {
            # Keep last state for a restart
            "s_track_upper_restart": (["lat", "lon"],s_track_upper),
            "s_track_lower_restart": (["lat", "lon"], s_track_lower),
            "s_track_upper": (["lat", "lon"], s_track_upper_mean),
            "s_track_lower": (["lat", "lon"], s_track_lower_mean),
            "e_track": (["lat", "lon"], e_track),
            "north_loss": (["lon"], north_loss),
            "south_loss": (["lon"], south_loss),
            "east_loss": (["lat"], east_loss),
            "west_loss": (["lat"],west_loss),
        }
    )
    return (s_track_upper, s_track_lower, ds)


region = xr.open_dataset(config["region"]).region_flood.values

for i, date in enumerate(reversed(datelist[:])):
    print(date)
    preprocessed_data = xr.open_dataset(input_path(date))

    # Resample to (higher) target frequency
    # After this, the fluxes will be "in between" the states
    fluxes, states = resample(preprocessed_data, config['target_frequency'])

    # Convert flux data to volumes
    change_units(fluxes, config["target_frequency"])

    # Apply a stability correction if needed
    stabilize_fluxes(fluxes, states)

    # Determine the vertical moisture flux
    fluxes["f_vert"] = calculate_fv(fluxes, states, config["periodic_boundary"], config["kvf"])

    if i == 0:
        if config["restart"]:
            # Reload last state from existing output
            ds = xr.open_dataset(output_path(date + pd.Timedelta(days=1)))
            s_track_upper = ds.s_track_upper_restart.values
            s_track_lower = ds.s_track_lower_restart.values
        else:
            # Allocate empty arrays based on shape of input data
            s_track_upper = np.zeros_like(states.s_upper[0])
            s_track_lower = np.zeros_like(states.s_upper[0])

    (s_track_upper, s_track_lower, processed_data) = backtrack(
        date,
        fluxes,
        states,
        s_track_upper,
        s_track_lower,
        region,
        config["kvf"],
    )

    # Write output to file
    # TODO: add (and cleanup) coordinates and units
    processed_data.to_netcdf(output_path(date))
