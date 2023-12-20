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
from wam2layers.tracking.io import load_data, load_region, output_path, write_output
from wam2layers.utils.profiling import ProgressTracker

logger = logging.getLogger(__name__)


def forwardtrack(
    F,
    S1,
    S0,
    region,
    output,
    config,
):
    # Unpack input data
    fx_upper = stagger_x(F["fx_upper"].values)
    fy_upper = stagger_y(F["fy_upper"].values)
    fx_lower = stagger_x(F["fx_lower"].values)
    fy_lower = stagger_y(F["fy_lower"].values)
    evap = F["evap"].values
    precip = F["precip"].values
    f_vert = F["f_vert"].values
    s_upper = S0["s_upper"].values
    s_lower = S0["s_lower"].values
    
    s_track_upper = output["s_track_upper_restart"].values
    s_track_lower = output["s_track_lower_restart"].values

    # Short name for often used expressions
    precip_lower = precip * s_lower / (s_upper + s_lower)
    precip_upper = precip * s_upper / (s_upper + s_lower)
    s_track_relative_lower = np.minimum(s_track_lower / s_lower, 1.0)
    s_track_relative_upper = np.minimum(s_track_upper / s_upper, 1.0)

    # Actual tracking (note: backtracking, all fluxes change sign)
    bc = config.periodic_boundary  # boundary condition True/False
    # TODO: apply terms in successive steps instead of all at once?
    s_track_lower += (
        horizontal_advection(s_track_relative_lower, fx_lower, fy_lower, bc)
        + vertical_advection(f_vert, s_track_relative_lower, s_track_relative_upper)
        + vertical_dispersion(
            f_vert, s_track_relative_lower, s_track_relative_upper, config.kvf
        )
        + region * evap 
        - precip_lower * s_track_relative_lower
    )

    s_track_upper += (
        horizontal_advection(s_track_relative_upper, fx_upper, fy_upper, bc)
        - vertical_advection(f_vert, s_track_relative_lower, s_track_relative_upper)
        - vertical_dispersion(
            f_vert, s_track_relative_lower, s_track_relative_upper, config.kvf
        )
        - precip_upper * s_track_relative_upper
    )

    # down and top: redistribute unaccounted water that is otherwise lost from the sytem
    # TODO build in logging for lost moisture
    lower_to_upper = np.maximum(0, s_track_lower - S1["s_lower"])
    upper_to_lower = np.maximum(0, s_track_upper - S1["s_upper"])
    s_track_lower = np.minimum(s_track_lower - lower_to_upper + upper_to_lower, s_lower)
    s_track_upper = np.minimum(s_track_upper - upper_to_lower + lower_to_upper, s_upper)
 
    # Update output fields
    output["p_track_lower"] += precip_lower * s_track_relative_lower
    output["p_track_upper"] += precip_upper * s_track_relative_upper
    
    # Bookkeep boundary losses as "tracked moisture at grid edges"
    output["p_track_lower"][0, :] += s_track_lower[0, :]
    output["p_track_lower"][-1, :] += s_track_lower[-1, :]
    output["p_track_upper"][0, :] += s_track_upper[0, :]
    output["p_track_upper"][-1, :] += s_track_upper[-1, :]
    s_track_upper[0, :] = 0
    s_track_upper[-1, :] = 0
    s_track_lower[0, :] = 0
    s_track_lower[-1, :] = 0
    if config.periodic_boundary is False: # bookkeep west and east losses
        output["p_track_lower"][:, 0] += s_track_lower[:, 0]
        output["p_track_lower"][:, -1] += s_track_lower[:, -1]
        output["p_track_upper"][:, 0] += s_track_upper[:, 0]
        output["p_track_upper"][:, -1] += s_track_upper[:, -1]
        s_track_upper[:, 0] = 0
        s_track_upper[:, -1] = 0
        s_track_lower[:, 0] = 0
        s_track_lower[:, -1] = 0

    output["s_track_lower_restart"].values = s_track_lower
    output["s_track_upper_restart"].values = s_track_upper
    output["tagged_evap"] += region * evap


def initialize(config_file):
    """Read config, region, and initial states."""
    logger.info(f"Initializing experiment with config file {config_file}")

    config = Config.from_yaml(config_file)
    region = load_region(config)

    output = initialize_outputs(region)
    region = region.values

    if config.restart:
        # Reload last state from existing output
        date = config.track_end_date
        ds = xr.open_dataset(output_path(date, config))
        output["s_track_upper_restart"].values = ds.s_track_upper_restart.values
        output["s_track_lower_restart"].values = ds.s_track_lower_restart.values

    logger.info(f"Output will be written to {config.output_folder.absolute()}.")
    return config, region, output


def initialize_outputs(region):
    """Allocate output arrays."""

    proto = region
    output = xr.Dataset(
        {
            # Keep last state for a restart
            "s_track_upper_restart": xr.zeros_like(proto),
            "s_track_lower_restart": xr.zeros_like(proto),
            "p_track_upper": xr.zeros_like(proto),
            "p_track_lower": xr.zeros_like(proto),
            "tagged_evap": xr.zeros_like(proto),
        }
    )
    return output


def initialize_time(config, direction="forward"):
    dt = pd.Timedelta(config.target_frequency)

    if direction == "forward":
        t0 = pd.Timestamp(config.track_start_date)
        th = t0 + dt / 2  # th is "half" time, i.e. between t0 and t1
        t1 = t0 + dt
    elif direction == "backward":
        t1 = pd.Timestamp(config.track_end_date)
        th = t1 - dt / 2  # th is "half" time, i.e. between t0 and t1
        t0 = t1 - dt
    else:
        raise ValueError("Direction should be forward or backward")

    return t0, th, t1, dt


#############################################################################
# With the correct import statements, the code in the function below could
# alternatively be be used as a script in a separate python file or notebook.
#############################################################################


def run_experiment(config_file):
    """Run a forwardtracking experiment from start to finish."""
    config, region, output = initialize(config_file)

    event_start, event_end = config.event_start_date, config.event_end_date
    progress_tracker = ProgressTracker(output, mode="forwardtrack")

    t0, th, t1, dt = initialize_time(config, direction="forward")

    while t1 <= config.track_end_date:
        S0 = load_data(t0, config, "states")
        F = load_data(th, config, "fluxes")
        S1 = load_data(t1, config, "states")

        # Convert data to volumes
        change_units(S0, config.target_frequency)
        change_units(F, config.target_frequency)
        change_units(S1, config.target_frequency)

        # Apply a stability correction if needed
        stabilize_fluxes(F, S0, progress_tracker, config, t0)

        # Determine the vertical moisture flux
        F["f_vert"] = calculate_fz(F, S0, S1, config.kvf)

        # Only track the precipitation at certain timesteps
        if not event_start <= t0 <= event_end:
            F["evap"] = 0

        # TODO: consider changing signature to F, S0, S1
        # Inside forwardtrack the "output" dictionary is updated
        forwardtrack(F, S1, S0, region, output, config)
        t1 += dt
        th += dt
        t0 += dt

        # Daily output
        is_output_time = t0 == t0.floor(config.output_frequency)
        is_final_step = t1 > config.track_end_date
        if is_output_time or is_final_step:
            progress_tracker.print_progress(t0, output)
            progress_tracker.store_intermediate_states(output)
            write_output(output, t0, config)

    logger.info("Experiment complete.")


###########################################################################
# The code below makes it possible to run wam2layers from the command line:
# >>> python forwardtrack.py path/to/cases/era5_2021.yaml
# or even:
# >>> wam2layers forwardtrack path/to/cases/era5_2021.yaml
###########################################################################


@click.command()
@click.argument("config_file", type=click.Path(exists=True))
def cli(config_file):
    """Run WAM2layers forwardtrack experiment from the command line.

    CONFIG_FILE: Path to WAM2layers experiment configuration file.

    Usage examples:

        \b
        - python path/to/forwardtrack.py path/to/cases/era5_2021.yaml
        - wam2layers forwardtrack path/to/cases/era5_2021.yaml
    """
    logger.info("Welcome to WAM2layers.")
    logger.info("Starting forwardtrack experiment.")
    run_experiment(config_file)


if __name__ == "__main__":
    cli()
