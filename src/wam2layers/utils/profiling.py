import logging
import time

import numpy as np
import pandas as pd
import psutil
import xarray as xr

from wam2layers.config import Config

logger = logging.getLogger(__name__)


class Profiler:
    def __init__(self):
        self.t0 = time.time()
        self.mem0 = self._current_memory()

    def __call__(self):
        """Report memory and time increments."""
        dt = round((time.time() - self.t0), 2)
        dmem = self._current_memory() - self.mem0
        # self.t0 = time.time()
        # self.mem0 = self._current_memory()
        return dt, dmem

    def _current_memory(self):
        return psutil.Process().memory_info().rss / (1024 * 1024)


class ProgressTracker:
    def __init__(self, output, mode="backtrack"):
        """Keep track of tagged and tracked moisture."""
        self.mode = mode
        self.total_tagged_moisture = 0
        self.tracked = 0
        self.lost_water = 0
        self.gained_water = 0
        self.boundary_transport = 0
        self.store_intermediate_states(output)
        self.profile = Profiler()
        self.stability_correction_previous_grid = 0.0
        self.stability_correction_previous_factor = 1.0

    def store_intermediate_states(self, output):
        """Update moisture totals based on output since previous time step.

        Intended for when the output fields are reset at write time, so we
        don't lose accumulations from previous outputs.
        """
        totals = output.sum()

        if self.mode == "backtrack":
            self.tracked += totals["e_track"]
            self.total_tagged_moisture += totals["tagged_precip"]
        else:  # mode is forwardtrack
            self.tracked += totals["p_track_lower"] + totals["p_track_upper"]
            self.total_tagged_moisture += totals["tagged_evap"]

        # Don't include boundary losses in log messages
        interior_losses = output.losses[1:-1, 1:-1].sum()
        boundary_transport = totals["losses"] - interior_losses

        self.lost_water += interior_losses
        self.boundary_transport += boundary_transport
        self.gained_water += totals["gains"]

    def print_progress(self, t, output):
        """Print some useful runtime diagnostics."""
        totals = output.sum()

        if self.mode == "backtrack":
            tracked = self.tracked + totals["e_track"]
            total_tagged_moisture = self.total_tagged_moisture + totals["tagged_precip"]
        else:  # mode is forwardtrack
            tracked = self.tracked + totals["p_track_upper"] + totals["p_track_lower"]
            total_tagged_moisture = self.total_tagged_moisture + totals["tagged_evap"]
        still_in_atmosphere = (
            totals["s_track_upper_restart"] + totals["s_track_lower_restart"]
        )

        # Don't include boundary losses in log messages
        interior_losses = output.losses[1:-1, 1:-1].sum()
        boundary_transport = totals["losses"] - interior_losses

        lost_water = self.lost_water + interior_losses
        boundary_transport = self.boundary_transport + boundary_transport
        gained_water = self.gained_water + totals["gains"]

        tracked_percentage = tracked / total_tagged_moisture * 100
        in_atmos_percentage = still_in_atmosphere / total_tagged_moisture * 100
        lost_percentage = lost_water / total_tagged_moisture * 100
        gained_percentage = gained_water / total_tagged_moisture * 100
        boundary_percentage = boundary_transport / total_tagged_moisture * 100
        closure_check_percentage = (
            tracked_percentage
            + in_atmos_percentage
            + lost_percentage
            + gained_percentage
            + boundary_percentage
        )

        time, memory = self.profile()

        # TODO: print boundary losses and internal losses separately
        logger.info(
            f"{t} - "
            f"Tracked moisture: {tracked_percentage.item():.2f}%. "
            f"Tagged in atmosphere: {in_atmos_percentage.item():.2f}%. "
            f"Boundary transport: {boundary_percentage.item():.2f}%. "
            f"Lost moisture: {lost_percentage.item():.2f}%. "
            f"Gained moisture: {gained_percentage.item():.2f}%. "
            f"Tagging closure: {closure_check_percentage.item():.2f}%. "
            f"Time since start: {time}s, RAM: {memory:.2f} MB"
        )

    def track_stability_correction(
        self,
        fy_corrected: xr.DataArray,
        fy_abs: xr.DataArray,
        config: Config,
        t: pd.Timestamp,
    ) -> None:
        """Issue warning if correction exceeds criterion.

        Warning advises to reduce the timestep.

        Criteria:
        1. A correction is applied for more than 5% of the grid (first time)
        2. At least 10% more of the grid is corrected compared to the last warning
        3. Previous strongest correction factor is exceeded by 10%
        """
        corrected_area = fy_corrected < fy_abs
        corrected_percent = corrected_area.sum() / corrected_area.count() * 100
        correction_factor = fy_corrected / fy_abs
        strongest_correction = correction_factor.min()

        warn = False
        if strongest_correction < self.stability_correction_previous_factor / 1.1:
            self.stability_correction_previous_factor = strongest_correction
            warn = True

        elif (
            corrected_percent > 5  # at least 5 percent of the grid should be corrected
            and corrected_percent > self.stability_correction_previous_grid * 1.1
        ):
            self.stability_correction_previous_grid = corrected_percent
            warn = True

        if warn:
            # Write correction field to output debug file
            # TODO: doesn't feel like the right place to do this.
            # a model object could improve the code structure.
            debug_dir = config.output_folder / "debug"
            debug_dir.mkdir(exist_ok=True)
            timestamp = t.strftime("%Y%m%d-%H%M%S")
            filename = debug_dir / f"stability_correction_{timestamp}.nc"

            correction_factor.name = "correction_factor"
            correction_factor.to_netcdf(filename)

            logger.warning(
                f"Warning: stability correction applied to {corrected_percent:.1f}% of "
                "grid.\n"
                f"    Average correction factor was {correction_factor.mean():.2f}, "
                f"the strongest correction factor was {strongest_correction:.2f}\n"
                f"    The correction factor field is written to {filename}.\n"
                "    The next warning will only be raised when at least 10% more cells "
                "are corrected, or the strongest correction factor is 10% stronger."
            )
