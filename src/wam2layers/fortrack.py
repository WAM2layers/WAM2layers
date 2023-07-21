from functools import lru_cache

import click
import numpy as np
import pandas as pd
import xarray as xr

from wam2layers.config import Config
from wam2layers.preprocessing.shared import get_grid_info
from wam2layers.utils.profiling import ProgressTracker

from wam2layers.tracking.backtrack import get_tracking_dates
from wam2layers.tracking.backtrack import input_path

def output_path(date, config):
    output_dir = config.output_folder
    return f"{output_dir}/fortrack_{date.strftime('%Y-%m-%dT%H:%M')}.nc"

# LRU Cache keeps the file open so we save a bit on I/O
@lru_cache(maxsize=2)
from wam2layers.tracking.backtrack import read_data_at_date

# This one can't be cached as we'll be overwriting the content.
from wam2layers.tracking.backtrack import read_data_at_time

from wam2layers.tracking.backtrack import load_data

def load_fluxes(t):
    t_current = t + pd.Timedelta(config.target_frequency) / 2
    return load_data(t_current, "fluxes")

def load_states(t):
    t_prev = t
    t_next = t + pd.Timedelta(config.target_frequency)
    states_prev = load_data(t_prev, "states")
    states_next = load_data(t_next, "states")
    return states_prev, states_next

from wam2layers.tracking.backtrack import time_in_range

from wam2layers.tracking.backtrack import load_region

from wam2layers.tracking.backtrack import to_edges_zonal
from wam2layers.tracking.backtrack import to_edges_meridional

from wam2layers.tracking.backtrack import look_north
from wam2layers.tracking.backtrack import look_south
from wam2layers.tracking.backtrack import look_east
from wam2layers.tracking.backtrack import look_west

from wam2layers.tracking.backtrack import split_vertical_flux
from wam2layers.tracking.backtrack import change_units

from wam2layers.tracking.backtrack import stabilize_fluxes

from wam2layers.tracking.backtrack import convergence

from wam2layers.tracking.backtrack import convergence

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

    # Determine horizontal fluxes over the grid-cell boundaries
    f_e_lower_we, f_e_lower_ew, f_w_lower_we, f_w_lower_ew = to_edges_zonal(
        fx_lower, config.periodic_boundary
    )
    f_e_upper_we, f_e_upper_ew, f_w_upper_we, f_w_upper_ew = to_edges_zonal(
        fx_upper, config.periodic_boundary
    )

    (
        fy_n_lower_sn,
        fy_n_lower_ns,
        fy_s_lower_sn,
        fy_s_lower_ns,
    ) = to_edges_meridional(fy_lower)
    (
        fy_n_upper_sn,
        fy_n_upper_ns,
        fy_s_upper_sn,
        fy_s_upper_ns,
    ) = to_edges_meridional(fy_upper)

    # Short name for often used expressions
    s_track_relative_lower = s_track_lower / s_lower
    s_track_relative_upper = s_track_upper / s_upper
    if config.periodic_boundary:
        inner = np.s_[1:-1, :]
    else:
        inner = np.s_[1:-1, 1:-1]

    # Actual tracking (note: backtracking, all terms have been negated)
    s_track_lower[inner] += (
        +f_e_lower_we * look_east(s_track_relative_lower)
        + f_w_lower_ew * look_west(s_track_relative_lower)
        + fy_n_lower_sn * look_north(s_track_relative_lower)
        + fy_s_lower_ns * look_south(s_track_relative_lower)
        + f_upward * s_track_relative_upper
        - f_downward * s_track_relative_lower
        - fy_s_lower_sn * s_track_relative_lower
        - fy_n_lower_ns * s_track_relative_lower
        - f_e_lower_ew * s_track_relative_lower
        - f_w_lower_we * s_track_relative_lower
        + tagged_precip * (s_lower / s_total)
        - evap * s_track_relative_lower
    )[inner]

    s_track_upper[inner] += (
        +f_e_upper_we * look_east(s_track_relative_upper)
        + f_w_upper_ew * look_west(s_track_relative_upper)
        + fy_n_upper_sn * look_north(s_track_relative_upper)
        + fy_s_upper_ns * look_south(s_track_relative_upper)
        + f_downward * s_track_relative_lower
        - f_upward * s_track_relative_upper
        - fy_s_upper_sn * s_track_relative_upper
        - fy_n_upper_ns * s_track_relative_upper
        - f_w_upper_we * s_track_relative_upper
        - f_e_upper_ew * s_track_relative_upper
        + tagged_precip * (s_upper / s_total)
    )[inner]

    # down and top: redistribute unaccounted water that is otherwise lost from the sytem
    lower_to_upper = np.maximum(0, s_track_lower - states_next["s_lower"])
    upper_to_lower = np.maximum(0, s_track_upper - states_next["s_upper"])
    s_track_lower[inner] = (s_track_lower - lower_to_upper + upper_to_lower)[inner]
    s_track_upper[inner] = (s_track_upper - upper_to_lower + lower_to_upper)[inner]

    # Update output fields
    output["e_track"] += evap * (s_track_lower / s_lower)
    output["north_loss"] += (
        fy_n_upper_ns * s_track_relative_upper + fy_n_lower_ns * s_track_relative_lower
    )[1, :]
    output["south_loss"] += (
        fy_s_upper_sn * s_track_relative_upper + fy_s_lower_sn * s_track_relative_lower
    )[-2, :]

    if config.periodic_boundary == False:
        output["east_loss"] += (
            f_e_upper_ew * s_track_relative_upper
            + f_e_lower_ew * s_track_relative_lower
        )[:, -2]
        output["west_loss"] += (
            f_w_upper_we * s_track_relative_upper
            + f_w_lower_we * s_track_relative_lower
        )[:, 1]

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
