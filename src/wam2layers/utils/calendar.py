"""Functionality related to (non-standard) calendars."""
import warnings
from pathlib import Path

import cftime
import xarray as xr

from wam2layers.config import Config


def has_cftime(ds: xr.Dataset) -> bool:
    """Check if the dataset has dimension 'time' with cftime values."""
    return "cftime" in str(type(ds["time"].values[0]))


def get_calendar_type(ds: xr.Dataset) -> str:
    """Look up the calendar type using cftime."""
    if has_cftime(ds):
        return ds["time"].values[0].calendar
    else:
        return "standard"


def validate_calendar_type(config: Config) -> None:
    file = template_to_files(config)[0]
    data_calendar = get_calendar_type(xr.open_dataset(file))
    if data_calendar != config.calendar:
        warnings.warn(
            "Calendar mismatch detected:\n"
            f"The configuration file has calendar type '{config.calendar}', "
            f"while the data has calendar '{data_calendar}'.\n"
            "This can lead to unpredictable errors!"
        )


def template_to_files(config: Config, var: str) -> list[Path]:
    template = Path(config.filename_template.format(variable=var))

    return list(template.parent.glob(template.name))
