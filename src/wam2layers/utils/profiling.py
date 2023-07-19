import time
import numpy as np

import psutil
import logging

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
    def __init__(self, output):
        """Keep track of tagged and tracked moisture."""
        self.total_tagged_moisture = 0
        self.tracked = 0
        self.boundary = 0
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
        self.tracked += totals["e_track"]
        self.total_tagged_moisture += totals["tagged_precip"]
        self.boundary += (
            totals["north_loss"]
            + totals["south_loss"]
            + totals["east_loss"]
            + totals["west_loss"]
        )

    def print_progress(self, t, output):
        """Print some useful runtime diagnostics."""
        totals = output.sum()
        tracked = self.tracked + totals["e_track"]
        total_tagged_moisture = self.total_tagged_moisture + totals["tagged_precip"]
        boundary = self.boundary + (
            totals["north_loss"]
            + totals["south_loss"]
            + totals["east_loss"]
            + totals["west_loss"]
        )
        still_in_atmosphere = (
            totals["s_track_upper_restart"] + totals["s_track_lower_restart"]
        )

        total_tracked_moisture = tracked + still_in_atmosphere + boundary
        tracked_percentage = (tracked + boundary) / total_tagged_moisture * 100
        lost_percentage = (1 - total_tracked_moisture / total_tagged_moisture) * 100

        time, memory = self.profile()

        logger.info(
            f"{t} - "
            f"Tracked moisture: {tracked_percentage.item():.2f}%. "
            f"Lost moisture: {lost_percentage.item():.2f}%. "
            f"Time since start: {time}s, RAM: {memory:.2f} MB"
        )

    def track_stability_correction(self, fy_corrected, fy_abs):
        """Issue warning if correction exceeds criterion.

        Warning advises to reduce the timestep.

        Criteria:
        1. A correction is applied for more than 5% of the grid (first time)
        2. Previous correction is exceeded by either 5%-point of grid
        3. Previous correction is doubled in magnitude
        """
        corrected = fy_corrected < fy_abs
        corrected_percent = corrected.sum() / corrected.count() * 100
        correction = np.where(corrected, fy_abs - fy_corrected, 0).mean()  # todo MAX?

        # Reversed conditions lead to cleaner code
        if correction < (2 * self.stability_correction_previous_value):
            return
        if (corrected_percent - 5) < self.stability_correction_previous_value:
            return

        self.stability_correction_previous_grid = corrected_percent
        self.stability_correction_previous_value = correction

        logger.warn(
            f"Stability correction applied to {corrected_percent:.1f}% of "
            f" grid, average correction was {correction.mean():.1f}"
        )
