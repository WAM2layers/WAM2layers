from functools import lru_cache
import functools

import click
import numpy as np
import pandas as pd
import xarray as xr

from wam2layers.config import Config
from wam2layers.preprocessing.shared import get_grid_info
from wam2layers.utils.profiling import ProgressTracker


def get_tracking_dates(config):
    """Dates for tracking."""
    # E.g. if data from 00:00 to 23:00
    # Using target freq directly would end at 23:45 or so
    input_dates = pd.date_range(
        start=config.track_start_date,
        end=config.track_end_date,
        freq=config.input_frequency,
    )

    return pd.date_range(
        start=input_dates[0],
        end=input_dates[-1],
        freq=config.target_frequency,
        inclusive="right",
    )


def input_path(date, config):
    input_dir = config.preprocessed_data_folder
    return f"{input_dir}/{date.strftime('%Y-%m-%d')}_fluxes_storages.nc"


def output_path(date, config):
    output_dir = config.output_folder
    return f"{output_dir}/backtrack_{date.strftime('%Y-%m-%dT%H-%M')}.nc"


# LRU Cache keeps the file open so we save a bit on I/O
@lru_cache(maxsize=2)
def read_data_at_date(d):
    """Load input data for given date."""
    file = input_path(d, config)
    return xr.open_dataset(file, cache=False)


# This one can't be cached as we'll be overwriting the content.
def read_data_at_time(t):
    """Get a single time slice from input data at time t."""
    ds = read_data_at_date(t)
    return ds.sel(time=t, drop=True)


def load_data(t, subset="fluxes"):
    """Load variable at t, interpolate if needed."""
    variables = {
        "fluxes": ["fx_upper", "fx_lower", "fy_upper", "fy_lower", "evap", "precip"],
        "states": ["s_upper", "s_lower"],
    }

    t1 = t.ceil(config.input_frequency)
    da1 = read_data_at_time(t1)[variables[subset]]
    if t == t1:
        # this saves a lot of work if times are aligned with input
        return da1

    t0 = t.floor(config.input_frequency)
    da0 = read_data_at_time(t0)[variables[subset]]
    if t == t0:
        return da0

    return da0 + (t - t0) / (t1 - t0) * (da1 - da0)


def load_fluxes(t):
    t_current = t - pd.Timedelta(config.target_frequency) / 2
    return load_data(t_current, "fluxes")


def load_states(t):
    t_prev = t
    t_next = t - pd.Timedelta(config.target_frequency)
    states_prev = load_data(t_prev, "states")
    states_next = load_data(t_next, "states")
    return states_prev, states_next


def time_in_range(start, end, current):
    """Returns whether current is in the range [start, end]"""
    return start <= current <= end


def load_region(config):
    # TODO: make variable name more generic
    return xr.open_dataset(config.region).source_region


def stagger_x(f):
    """Interpolate f to the grid cell interfaces.

    Only the values at the interior interfaces are returned

    Arguments:
        f: 2d array of shape [M, N]

    Returns:
        2d array of shape [M-2, N-1]
    """
    return  0.5 * (f[:, :-1] + f[:, 1:])[1:-1, :]


def stagger_y(f):
    """Interpolate f to the grid cell interfaces.

    Only the values at the interior interfaces are returned

    Arguments:
        f: 2d array of shape [M, N]

    Returns:
        2d array of shape [M-1, N-2]
    """
    return 0.5 * (f[:-1, :] + f[1:, :])[:, 1:-1]


def pad_boundaries(*args, periodic=False):
    """Add boundary padding to all input arrays."""
    periodic_x = functools.partial(np.pad, pad_width=((0, 0), (1, 1)), mode='wrap')
    zero_y = functools.partial(np.pad, pad_width=((1, 1), (0, 0)), mode='constant', constant_values=0)
    zero_xy = functools.partial(np.pad, pad_width=1, mode='constant', constant_values=0)
    if periodic:
        return [periodic_x(zero_y(arg)) for arg in args]
    return [zero_xy(arg) for arg in args]


def advection(q, u, v):
        """Calculate advection on a staggered grid.

        Can only calculate advection for the interior of the array. Hence the
        resulting array is 2 cells smaller in both directions.

        Arguments:
            q: array of shape [M, N]
            u: array of shape = [M-2, N-1]
            v: array of shape = [M-1, N-2]

        Returns:
            array of shape [M-2, N-2]

        Example:

            >>> q = np.zeros((5, 5))
            >>> q[2, 2] = 1
            >>> u = np.ones((3, 4)) * .5
            >>> v = np.ones((4, 3)) * .5
            >>> advection(q, u, v)
            array([[ 0. ,  0.5,  0. ],
                   [ 0. , -1. ,  0.5],
                   [ 0. ,  0. ,  0. ]])

        """
        west = np.s_[:-1]
        east = np.s_[1:]

        # TODO revert direction of era5 latitude in pre-processing(?)
        south = np.s_[1:]
        north = np.s_[:-1]

        inner = np.s_[1:-1]

        # Donor cell upwind scheme (2 directions seperately)
        uq = np.where(
            u > 0,
            u * q[inner, west],
            u * q[inner, east]
            )  # [M-2, N-1]

        vq = np.where(
            v > 0,
            v * q[south, inner],
            v * q[north, inner]
        )  # [M-1, N-2]

        adv_x = uq[:, west] - uq[:, east]  # [M-2, N-2]
        adv_y = vq[south, :] - vq[north, :]  # [M-2, N-2]

        return adv_x + adv_y  # [M-2, N-2]


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


def change_units(data, target_freq):
    """Change units to m3.
    Multiply by edge length or area to get flux in m3
    Multiply by time to get accumulation instead of flux
    Divide by density of water to go from kg to m3
    """
    density = 1000  # [kg/m3]
    a, ly, lx = get_grid_info(data)

    total_seconds = pd.Timedelta(target_freq).total_seconds()

    for variable in data.data_vars:
        if variable in ["fx_upper", "fx_lower"]:
            data[variable] *= total_seconds / density * ly
        elif variable in ["fy_upper", "fy_lower"]:
            data[variable] *= total_seconds / density * lx[:, None]
        elif variable in ["evap", "precip"]:
            data[variable] *= total_seconds / density * a[:, None]
        elif variable in ["s_upper", "s_lower"]:
            data[variable] *= a[:, None] / density
        else:
            raise ValueError(f"Unrecognized variable {variable}")
        data[variable] = data[variable].assign_attrs(units="m**3")


def stabilize_fluxes(current, previous):
    """Stabilize the outfluxes / influxes.

    CFL: Water cannot move further than one grid cell per timestep.
    """
    for level in ["upper", "lower"]:
        fx = current["fx_" + level]
        fy = current["fy_" + level]
        s = previous["s_" + level]

        fx_abs = np.abs(fx)
        fy_abs = np.abs(fy)
        ft_abs = fx_abs + fy_abs

        fx_corrected = 1 / 2 * fx_abs / ft_abs * s.values
        fx_stable = np.minimum(fx_abs, fx_corrected)

        fy_corrected = 1 / 2 * fy_abs / ft_abs * s.values
        fy_stable = np.minimum(fy_abs, fy_corrected)

        # Get rid of any nan values
        fx_stable.fillna(0)
        fy_stable.fillna(0)

        # Re-instate the sign
        current["fx_" + level] = np.sign(fx) * fx_stable
        current["fy_" + level] = np.sign(fy) * fy_stable


def convergence(fx, fy):
    # Note: latitude decreasing, hence positive fy gradient is convergence
    return np.gradient(fy, axis=-2) - np.gradient(fx, axis=-1)


def calculate_fv(fluxes, states_prev, states_next):
    """Calculate the vertical fluxes.

    Note: fluxes are given at temporal midpoints between states.
    """
    s_diff = states_prev - states_next
    s_mean = (states_prev + states_next) / 2
    s_total = s_mean.s_upper + s_mean.s_lower
    s_rel = s_mean / s_total

    tendency_upper = (
        convergence(fluxes.fx_upper, fluxes.fy_upper)
        - fluxes.precip.values * s_rel.s_upper
    )
    tendency_lower = (
        convergence(fluxes.fx_lower, fluxes.fy_lower)
        - fluxes.precip.values * s_rel.s_lower
        + fluxes.evap
    )

    residual_upper = s_diff.s_upper - tendency_upper
    residual_lower = s_diff.s_lower - tendency_lower

    # compute the resulting vertical moisture flux; the vertical velocity so
    # that the new residual_lower/s_lower = residual_upper/s_upper (positive downward)
    fv = s_rel.s_lower * (residual_upper + residual_lower) - residual_lower

    # stabilize the outfluxes / influxes; during the reduced timestep the
    # vertical flux can maximally empty/fill 1/x of the top or down storage
    stab = 1.0 / (config.kvf + 1.0)
    flux_limit = np.minimum(s_mean.s_upper, s_mean.s_lower)
    fv_stable = np.minimum(np.abs(fv), stab * flux_limit)

    # Reinstate the sign
    return np.sign(fv) * fv_stable


def backtrack(
    fluxes,
    states_prev,
    states_next,
    region,
    output,
):
    # Unpack input data
    fx_upper = fluxes["fx_upper"].values
    fy_upper = fluxes["fy_upper"].values
    fx_lower = fluxes["fx_lower"].values
    fy_lower = fluxes["fy_lower"].values
    evap = fluxes["evap"].values
    precip = fluxes["precip"].values
    f_vert = fluxes["f_vert"].values
    s_upper = states_prev["s_upper"].values
    s_lower = states_prev["s_lower"].values

    s_track_upper = output["s_track_upper_restart"].values
    s_track_lower = output["s_track_lower_restart"].values

    tagged_precip = region * precip
    s_total = s_upper + s_lower

    # separate the direction of the vertical flux and make it absolute
    f_downward, f_upward = split_vertical_flux(config.kvf, f_vert)

    # Determine fluxes at the faces of the grid cells
    fx_lower = stagger_x(fx_lower)
    fx_upper = stagger_x(fx_upper)
    fy_lower = stagger_y(fy_lower)
    fy_upper = stagger_y(fy_upper)

    # TODO move staggering to preprocessing (?)

    # Short name for often used expressions
    s_track_relative_lower = np.minimum(s_track_lower / s_lower, 1.0)
    s_track_relative_upper = np.minimum(s_track_upper / s_upper, 1.0)

    # Actual tracking (note: backtracking, all terms have been negated)
    padded = functools.partial(pad_boundaries, periodic=config.periodic_boundary)
    s_track_lower += (
        + advection(*padded(s_track_relative_lower, -fx_lower, -fy_lower))
        + (f_upward * s_track_relative_upper)
        - (f_downward * s_track_relative_lower)
        + (tagged_precip * (s_lower / s_total))
        - (evap * s_track_relative_lower)
    )

    s_track_upper += (
        + advection(*padded(s_track_relative_upper, -fx_upper, -fy_upper))
        + (f_downward * s_track_relative_lower)
        - (f_upward * s_track_relative_upper)
        + (tagged_precip * (s_upper / s_total))
    )

    # down and top: redistribute unaccounted water that is otherwise lost from the sytem
    lower_to_upper = np.maximum(0, s_track_lower - states_next["s_lower"])
    upper_to_lower = np.maximum(0, s_track_upper - states_next["s_upper"])
    s_track_lower = (s_track_lower - lower_to_upper + upper_to_lower)
    s_track_upper = (s_track_upper - upper_to_lower + lower_to_upper)

    # Update output fields
    output["e_track"] += evap * np.minimum(s_track_lower / s_lower, 1.0)

    # Bookkeep boundary losses as "tracked moisture at grid edges"
    output["e_track"][0, :] += (s_track_upper + s_track_lower)[0, :]
    output["e_track"][-1, :] += (s_track_upper + s_track_lower)[-1, :]
    s_track_upper[0, :] = 0
    s_track_upper[-1, :] = 0
    s_track_lower[0, :] = 0
    s_track_lower[-1, :] = 0
    if config.periodic_boundary is False:
        output["e_track"][:, 0] += (s_track_upper + s_track_lower)[:, 0]
        output["e_track"][:, -1] += (s_track_upper + s_track_lower)[:, -1]
        s_track_upper[:, 0] = 0
        s_track_upper[:, -1] = 0
        s_track_lower[:, 0] = 0
        s_track_lower[:, -1] = 0

    output["s_track_lower_restart"].values = s_track_lower
    output["s_track_upper_restart"].values = s_track_upper
    output["tagged_precip"] += tagged_precip


def initialize(config_file):
    """Read config, region, and initial states."""
    print(f"Initializing experiment with config file {config_file}")

    config = Config.from_yaml(config_file)
    region = load_region(config)

    output = initialize_outputs(region)
    region = region.values

    if config.restart:
        # Reload last state from existing output
        tracking_dates = get_tracking_dates(config)
        date = tracking_dates[-1] + pd.Timedelta(days=1)
        ds = xr.open_dataset(output_path(date, config))
        output["s_track_upper_restart"].values = ds.s_track_upper_restart.values
        output["s_track_lower_restart"].values = ds.s_track_lower_restart.values

    print(f"Output will be written to {config.output_folder}.")
    return config, region, output


def initialize_outputs(region):
    """Allocate output arrays."""

    proto = region
    output = xr.Dataset(
        {
            # Keep last state for a restart
            "s_track_upper_restart": xr.zeros_like(proto),
            "s_track_lower_restart": xr.zeros_like(proto),
            "e_track": xr.zeros_like(proto),
            "tagged_precip": xr.zeros_like(proto),
            "north_loss": xr.zeros_like(proto.isel(latitude=0, drop=True)),
            "south_loss": xr.zeros_like(proto.isel(latitude=0, drop=True)),
            "east_loss": xr.zeros_like(proto.isel(longitude=0, drop=True)),
            "west_loss": xr.zeros_like(proto.isel(longitude=0, drop=True)),
        }
    )
    return output


def write_output(output, t):
    # TODO: add back (and cleanup) coordinates and units
    path = output_path(t, config)
    print(f"{t} - Writing output to file {path}")
    output.astype("float32").to_netcdf(path)

    # Flush previous outputs
    output[
        [
            "e_track",
            "tagged_precip",
            "north_loss",
            "south_loss",
            "east_loss",
            "west_loss",
        ]
    ] *= 0


#############################################################################
# With the correct import statements, the code in the function below could
# alternatively be be used as a script in a separate python file or notebook.
#############################################################################


def run_experiment(config_file):
    """Run a backtracking experiment from start to finish."""
    global config

    config, region, output = initialize(config_file)

    tracking_dates = get_tracking_dates(config)

    progress_tracker = ProgressTracker(output)
    for t in reversed(tracking_dates):
        fluxes = load_fluxes(t)
        states_prev, states_next = load_states(t)

        # Convert data to volumes
        change_units(states_prev, config.target_frequency)
        change_units(states_next, config.target_frequency)
        change_units(fluxes, config.target_frequency)

        # Apply a stability correction if needed
        stabilize_fluxes(fluxes, states_next)

        # Determine the vertical moisture flux
        fluxes["f_vert"] = calculate_fv(fluxes, states_prev, states_next)

        # Only track the precipitation at certain timesteps
        if (
            time_in_range(
                config.event_start_date,
                config.event_end_date,
                t,
            )
            == False
        ):
            fluxes["precip"] = 0

        backtrack(
            fluxes,
            states_prev,
            states_next,
            region,
            output,
        )

        # Daily output
        if t == t.floor(config.output_frequency) or t == tracking_dates[0]:
            progress_tracker.print_progress(t, output)
            progress_tracker.store_intermediate_states(output)
            write_output(output, t)

    print("Experiment complete.")


###########################################################################
# The code below makes it possible to run wam2layers from the command line:
# >>> python backtrack.py path/to/cases/era5_2021.yaml
# or even:
# >>> wam2layers backtrack path/to/cases/era5_2021.yaml
###########################################################################


@click.command()
@click.argument("config_file", type=click.Path(exists=True))
def cli(config_file):
    """Run WAM2layers backtrack experiment from the command line.

    CONFIG_FILE: Path to WAM2layers experiment configuration file.

    Usage examples:

        \b
        - python path/to/backtrack.py path/to/cases/era5_2021.yaml
        - wam2layers backtrack path/to/cases/era5_2021.yaml
    """
    print("Welcome to WAM2layers.")
    print("Starting backtrack experiment.")
    run_experiment(config_file)


if __name__ == "__main__":
    cli()
