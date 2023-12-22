import logging
import time
from datetime import datetime

import numpy as np
import psutil
import xarray as xr

from wam2layers.tracking.io import load_tagging_region

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
        self.store_intermediate_states(output)
        self.profile = Profiler()
        self.stability_correction_previous_grid = 0
        self.stability_correction_previous_value = 0

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
        self.gained_water += totals["gains"]
        self.lost_water += totals["losses"]

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
        lost_water = self.lost_water + totals["losses"]
        gained_water = self.gained_water + totals["gains"]

        tracked_percentage = tracked / total_tagged_moisture * 100
        in_atmos_percentage = still_in_atmosphere / total_tagged_moisture * 100
        lost_percentage = lost_water / total_tagged_moisture * 100
        gained_percentage = gained_water / total_tagged_moisture * 100
        closure_check_percentage = (
            tracked_percentage
            + in_atmos_percentage
            + lost_percentage
            - gained_percentage
        )

        time, memory = self.profile()

        # TODO: print boundary losses and internal losses separately
        logger.info(
            f"{t} - "
            f"Tracked moisture: {tracked_percentage.item():.2f}%. "
            f"Tagged in atmosphere: {in_atmos_percentage.item():.2f}%. "
            f"Lost moisture: {lost_percentage.item():.2f}%. "
            f"Gained moisture: {gained_percentage.item():.2f}%. "
            f"Tagging closure: {closure_check_percentage.item():.2f}%. "
            f"Time since start: {time}s, RAM: {memory:.2f} MB"
        )

    def track_stability_correction(self, fy_corrected, fy_abs, config, coords, t):
        """Issue warning if correction exceeds criterion.

        Warning advises to reduce the timestep.

        Criteria:
        1. A correction is applied for more than 5% of the grid (first time)
        2. Previous correction is exceeded by either 5%-point of grid
        3. Previous correction is doubled in magnitude
        """
        corrected = fy_corrected < fy_abs
        corrected_percent = corrected.sum() / corrected.count() * 100
        correction = np.where(corrected, fy_abs - fy_corrected, 0)
        correction_max = correction.max()

        # Reversed conditions lead to cleaner code
        if correction_max < (2 * self.stability_correction_previous_value):
            return
        if (corrected_percent - 5) < self.stability_correction_previous_value:
            return

        self.stability_correction_previous_grid = corrected_percent
        self.stability_correction_previous_value = correction_max

        # Write correction field to output debug file
        # TODO: doesn't feel like the right place to do this.
        # a model object could improve the code structure.
        debug_dir = config.output_folder / "debug"
        debug_dir.mkdir(exist_ok=True)
        timestamp = t.strftime("%Y%m%d-%H%M%S")
        filename = debug_dir / f"stability_correction_{timestamp}.nc"

        ncfile = xr.DataArray(correction, coords=coords, name="correction")
        ncfile.to_netcdf(filename)

        logger.warn(
            f"Stability correction applied to {corrected_percent:.1f}% of "
            f"grid. Average correction was {correction.mean():.1f}, "
            f"maximum correction was {correction.max()}. The total "
            f"correction field is written to {filename}."
        )
