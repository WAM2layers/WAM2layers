import logging
from pathlib import Path

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
from wam2layers.tracking.io import (
    load_data,
    load_region,
    load_tagging_region,
    output_path,
    write_output,
)
from wam2layers.utils.profiling import ProgressTracker

logger = logging.getLogger(__name__)


def backtrack(
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
    s_upper = S1["s_upper"].values
    s_lower = S1["s_lower"].values

    s_track_upper = output["s_track_upper_restart"].values
    s_track_lower = output["s_track_lower_restart"].values

    # Short name for often used expressions
    s_track_relative_lower = np.minimum(s_track_lower / s_lower, 1.0)
    s_track_relative_upper = np.minimum(s_track_upper / s_upper, 1.0)

    # Actual tracking (note: backtracking, all fluxes change sign)
    bc = config.periodic_boundary  # boundary condition True/False
    # TODO: apply terms in successive steps instead of all at once?
    s_track_lower += (
        +horizontal_advection(s_track_relative_lower, -fx_lower, -fy_lower, bc)
        + vertical_advection(-f_vert, s_track_relative_lower, s_track_relative_upper)
        + vertical_dispersion(
            -f_vert, s_track_relative_lower, s_track_relative_upper, config.kvf
        )
        + region * precip * s_lower / (s_upper + s_lower)
        - evap * s_track_relative_lower
    )

    s_track_upper += (
        +horizontal_advection(s_track_relative_upper, -fx_upper, -fy_upper, bc)
        - vertical_advection(-f_vert, s_track_relative_lower, s_track_relative_upper)
        - vertical_dispersion(
            -f_vert, s_track_relative_lower, s_track_relative_upper, config.kvf
        )
        + region * precip * s_upper / (s_upper + s_lower)
    )

    # down and top: redistribute unaccounted water that is otherwise lost from the sytem
    # TODO build in logging for lost moisture
    lower_to_upper = np.maximum(0, s_track_lower - S0["s_lower"])
    upper_to_lower = np.maximum(0, s_track_upper - S0["s_upper"])
    s_track_lower = np.minimum(
        s_track_lower - lower_to_upper + upper_to_lower, S0["s_lower"]
    )
    s_track_upper = np.minimum(
        s_track_upper - upper_to_lower + lower_to_upper, S0["s_upper"]
    )

    # Update output fields
    output["e_track"] += evap * s_track_relative_lower

    # Bookkeep boundary losses as "tracked moisture at grid edges"
    output["e_track"][0, :] += (s_track_upper + s_track_lower)[0, :]
    output["e_track"][-1, :] += (s_track_upper + s_track_lower)[-1, :]
    s_track_upper[0, :] = 0
    s_track_upper[-1, :] = 0
    s_track_lower[0, :] = 0
    s_track_lower[-1, :] = 0
    if config.periodic_boundary is False:  # bookkeep west and east losses
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

    # Initialize outputs as empty fields based on the input coords
    sample_ds = load_data(config.tracking_end_date, config, "states")
    field_coords = sample_ds.isel(time=0, drop=True).coords
    empty_field = xr.DataArray(0, coords=field_coords)
    output = xr.Dataset(
        {
            # Keep last state for a restart
            "s_track_upper_restart": empty_field,
            "s_track_lower_restart": empty_field,
            "e_track": empty_field,
            "tagged_precip": empty_field,
        }
    )

    if config.restart:
        # Reload last state from existing output
        date = config.tracking_end_date
        ds = xr.open_dataset(output_path(date, config, mode="backtrack"))
        output["s_track_upper_restart"].values = ds.s_track_upper_restart.values
        output["s_track_lower_restart"].values = ds.s_track_lower_restart.values

    logger.info(f"Output will be written to {config.output_folder.absolute()}.")
    return config, output


def initialize_time(config, direction="forward"):
    dt = pd.Timedelta(config.target_frequency)

    if direction == "forward":
        t0 = pd.Timestamp(config.tracking_start_date)
        th = t0 + dt / 2  # th is "half" time, i.e. between t0 and t1
        t1 = t0 + dt
    elif direction == "backward":
        t1 = pd.Timestamp(config.tracking_end_date)
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
    """Run a backtracking experiment from start to finish."""
    config, output = initialize(config_file)

    progress_tracker = ProgressTracker(output, mode="backtrack")

    t0, th, t1, dt = initialize_time(config, direction="backward")

    while t0 >= config.tracking_start_date:
        S0 = load_data(t0, config, "states")
        F = load_data(th, config, "fluxes")
        S1 = load_data(t1, config, "states")

        if isinstance(config.tagging_region, Path):
            region = load_tagging_region(config, t=t0).values
        else:
            # TODO make it clearer.
            bbox = config.tagging_region
            region = xr.where(
                cond=(
                    (bbox.south <= S0.latitude <= bbox.north)
                    & (
                        (bbox.west <= S0.longitude <= bbox.east)
                        | (
                            (bbox.west > bbox.east)
                            & (
                                (bbox.west <= S0.longitude)
                                | (S0.longitude <= bbox.east)
                            )
                        )
                    )
                ),
                x=1,
                y=0,
            ).values

        # Only tag the precipitation at certain timesteps
        if not config.tagging_start_date <= t1 <= config.tagging_end_date:
            # TODO don't do this if region is 3d
            region.fill(0)

        # Convert data to volumes
        change_units(S0, config.target_frequency)
        change_units(F, config.target_frequency)
        change_units(S1, config.target_frequency)

        # Apply a stability correction if needed
        stabilize_fluxes(F, S1, progress_tracker, config, t1)

        # Determine the vertical moisture flux
        F["f_vert"] = calculate_fz(F, S0, S1, config.kvf)

        # Inside backtrack the "output" dictionary is updated
        backtrack(F, S1, S0, region, output, config)
        t1 -= dt
        th -= dt
        t0 -= dt

        # Daily output
        is_output_time = t1 == t1.floor(config.output_frequency)
        is_final_step = t0 < config.tracking_start_date
        if is_output_time or is_final_step:
            progress_tracker.print_progress(t1, output)
            progress_tracker.store_intermediate_states(output)
            write_output(output, t1, config, mode="backtrack")
            # Flush previous outputs
            output[["e_track", "tagged_precip"]] *= 0

    logger.info("Experiment complete.")
