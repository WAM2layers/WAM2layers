import warnings
from pathlib import Path

import numpy as np
import yaml


def _warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    """A more compact warning format.

    https://stackoverflow.com/a/26433913
    """
    short_path = Path(filename).relative_to(Path(filename).parents[2])
    return f"{short_path}:{lineno}: {category.__name__}: {message}\n"


# Override the format with which warnings are printed
warnings.formatwarning = _warning_on_one_line


def check_monotonic_increase(data, dim="level"):
    """Check if data is monotonically increasing in the given dimension."""
    if not data.diff(dim).min() >= 0:
        warnings.warn(
            f"Data is not monotonically increasing in the dimension {dim} for {data.name}."
        )


def check_monotonic_decrease(data, dim="level"):
    """Check if data is monotonically decreasing in the given dimension."""
    if not data.diff(dim).max() <= 0:
        warnings.warn(
            f"Data is not monotonically decreasing in the dimension {dim} for {data.name}."
        )


def check_units(data, valid_units):
    """Check if data has the given units."""
    if "units" not in data.attrs:
        warnings.warn(f"{data.name} has no units")
    elif data.units not in valid_units:
        warnings.warn(
            f"{data.name} has units '{data.units}'; expected (one of) {valid_units}."
        )


def check_positive(data):
    """Check if data is positive."""
    if not data.min() >= 0:
        warnings.warn(f"Data for {data.name} is not positive")


def check_negative(data):
    """Check if data is negative."""
    if not data.max() <= 0:
        warnings.warn(f"Data for {data.name} is not negative.")


def check_nonnegative(data):
    """Check if data is nonnegative."""
    if not data.min() >= 0:
        warnings.warn(f"Data for {data.name} is not nonnegative.")


def check_order_of_magnitude(data, order):
    """Check if data is within the given order of magnitude."""
    if not data.min() >= 10 ** order and data.max() <= 10**order:
        warnings.warn(
            f"Data  for {data.name}is not within the order of magnitude {order}."
        )


def check_frequency(data, frequency):
    """Check if data has the given frequency."""
    if not data.time.dt.frequency == frequency:
        warnings.warn(
            f"Data  for {data.name} has frequency {data.time.dt.frequency} instead of {frequency}."
        )


def check_shape(data, ndim):
    if data.ndim != ndim:
        warnings.warn("Expected 3 dimensions, found {data.ndim} for {data.name}.")


def check_range(data, range):
    if not data.min() >= range[0] and data.max() <= range[1]:
        warnings.warn(f"Data for {data.name} is not within the range {min} to {max}.")


def check_uniform(coord):
    spacing = np.diff(coord)
    if not spacing.min() == spacing.max():
        warnings.warn(
            f"Coordinate spacing is not uniform for coord {coord} of {data.name}"
        )


def check_coords(data, coords):
    check_shape(data, ndim=len(coords))
    for coord in coords:
        if coord not in data.coords:
            warnings.warn(f"Data for {data.name} does not have a coordinate {coord}.")
        else:
            check_uniform(data[coord])
        if coord == "latitude":
            check_monotonic_decrease(data[coord], dim=coord)
        else:
            check_monotonic_increase(data[coord], dim=coord)


def check_input(data):
    for variable, var_info in VARIABLES.items():
        check_units(data[variable], var_info["valid_units"])
        check_coords(data[variable], var_info["coordinates"])
        check_range(data[variable], var_info["range"])
        if var_info["all_positive"]:
            check_positive(data[variable])


VARIABLES = yaml.safe_load(
    """
    evap:
        valid_units: ["kg m-2 s-1", "kg/m2/s"]
        range: [0, 0.001]
        all_positive: True
        coordinates: [time, latitude, longitude]
    precip:
        valid_units: ["kg m-2 s-1", "kg/m2/s"]
        range: [0, 0.1]
        all_positive: True
        coordinates: [time, latitude, longitude]
    s_upper:
        valid_units: ["kg m-2", "kg/m2"]
        range: [0, 1000]
        all_positive: True
        coordinates: [time, latitude, longitude]
    s_lower:
        valid_units: ["kg m-2", "kg/m2"]
        range: [0, 1000]
        all_positive: True
        coordinates: [time, latitude, longitude]
    fx_upper:
        valid_units: ["kg m-1 s-1", "kg m^-1 s^-1", "kg m**-1 s**-1"]
        range: [-1000, 1000]
        all_positive: False
        coordinates: [time, latitude, longitude]
    fy_upper:
        valid_units: ["kg m-1 s-1", "kg m^-1 s^-1", "kg m**-1 s**-1"]
        range: [-1000, 1000]
        all_positive: False
        coordinates: [time, latitude, longitude]
    fx_lower:
        valid_units: ["kg m-1 s-1", "kg m^-1 s^-1", "kg m**-1 s**-1"]
        range: [-1000, 1000]
        all_positive: False
        coordinates: [time, latitude, longitude]
    fy_lower:
        valid_units: ["kg m-1 s-1", "kg m^-1 s^-1", "kg m**-1 s**-1"]
        range: [-1000, 1000]
        all_positive: False
        coordinates: [time, latitude, longitude]
"""
)
