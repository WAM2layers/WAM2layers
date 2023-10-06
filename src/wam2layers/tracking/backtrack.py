import logging

import click
import numpy as np
import pandas as pd
import xarray as xr

from wam2layers.config import Config
from wam2layers.preprocessing.shared import (
    calculate_fz,
    change_units,
    stabilize_fluxes,
    stagger_x,
    stagger_y,
)
from wam2layers.tracking.core import (
    horizontal_advection,
    vertical_advection,
    vertical_dispersion,
)
from wam2layers.tracking.io import load_data, load_region, output_path
from wam2layers.utils.profiling import ProgressTracker

logger = logging.getLogger(__name__)


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


def backtrack(
    F,
    S1,
    S0,
    region,
    output,
):
    # Unpack input data
    # TODO move staggering to preprocessing (?)
    fx_upper = stagger_x(F["fx_upper"].values)
    fy_upper = stagger_y(F["fy_upper"].values)
    fx_lower = stagger_x(F["fx_lower"].values)
    fy_lower = stagger_y(F["fy_lower"].values)
    evap = F["evap"].values
    precip = F["precip"].values
    f_vert = F["f_vert"].values
    s_upper = S1["s_upper"].values
    s_lower = S1["s_lower"].values

    s_track_upper = output["s_track_upper_restart"].values
    s_track_lower = output["s_track_lower_restart"].values

    # Short name for often used expressions
    s_track_relative_lower = np.minimum(s_track_lower / s_lower, 1.0)
    s_track_relative_upper = np.minimum(s_track_upper / s_upper, 1.0)

    # Actual tracking (note: backtracking, all fluxes change sign)
    bc = config.periodic_boundary
    s_track_lower += (
        +horizontal_advection(s_track_relative_lower, -fx_lower, -fy_lower, bc)
        + vertical_advection(-f_vert, s_track_relative_lower, s_track_relative_upper)
        + vertical_dispersion(-f_vert, s_track_relative_lower, s_track_relative_upper)
        + region * precip * s_lower / (s_upper + s_lower)
        - evap * s_track_relative_lower
    )

    s_track_upper += (
        +horizontal_advection(s_track_relative_upper, -fx_upper, -fy_upper, bc)
        - vertical_advection(-f_vert, s_track_relative_lower, s_track_relative_upper)
        - vertical_dispersion(-f_vert, s_track_relative_lower, s_track_relative_upper)
        + region * precip * s_upper / (s_upper + s_lower)
    )

    # down and top: redistribute unaccounted water that is otherwise lost from the sytem
    lower_to_upper = np.maximum(0, s_track_lower - S0["s_lower"])
    upper_to_lower = np.maximum(0, s_track_upper - S0["s_upper"])
    s_track_lower = s_track_lower - lower_to_upper + upper_to_lower
    s_track_upper = s_track_upper - upper_to_lower + lower_to_upper

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
    output["tagged_precip"] += region * precip


def initialize(config_file):
    """Read config, region, and initial states."""
    logger.info(f"Initializing experiment with config file {config_file}")

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

    logger.info(f"Output will be written to {config.output_folder}.")
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
        }
    )
    return output


def write_output(output, t):
    # TODO: add back (and cleanup) coordinates and units
    path = output_path(t, config)
    logger.info(f"{t} - Writing output to file {path}")
    output.astype("float32").to_netcdf(path)

    # Flush previous outputs
    output[["e_track", "tagged_precip"]] *= 0


#############################################################################
# With the correct import statements, the code in the function below could
# alternatively be be used as a script in a separate python file or notebook.
#############################################################################


def run_experiment(config_file):
    """Run a backtracking experiment from start to finish."""
    global config

    config, region, output = initialize(config_file)

    tracking_dates = get_tracking_dates(config)
    event_start, event_end = config.event_start_date, config.event_end_date

    progress_tracker = ProgressTracker(output)
    for t1 in reversed(tracking_dates):
        t0 = t1 - pd.Timedelta(config.target_frequency)
        th = t1 - pd.Timedelta(config.target_frequency) / 2

        S1 = load_data(t1, "states")
        S0 = load_data(t0, "states")
        F = load_data(th, "fluxes")

        # Convert data to volumes
        change_units(S1, config.target_frequency)
        change_units(S0, config.target_frequency)
        change_units(F, config.target_frequency)

        # Apply a stability correction if needed
        stabilize_fluxes(F, S0, progress_tracker, t1)

        # Determine the vertical moisture flux
        F["f_vert"] = calculate_fz(F, S0, S1)

        # Only track the precipitation at certain timesteps
        if not event_start <= t1 <= event_end:
            F["precip"] = 0

        backtrack(F, S1, S0, region, output)

        # Daily output
        if t1 == t1.floor(config.output_frequency) or t1 == tracking_dates[0]:
            progress_tracker.print_progress(t1, output)
            progress_tracker.store_intermediate_states(output)
            write_output(output, t1)

    logger.info("Experiment complete.")


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
    logger.info("Welcome to WAM2layers.")
    logger.info("Starting backtrack experiment.")
    run_experiment(config_file)


if __name__ == "__main__":
    cli()
