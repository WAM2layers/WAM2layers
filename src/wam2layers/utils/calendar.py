"""Functionality related to (non-standard) calendars."""

import re
import warnings
from pathlib import Path
from typing import TYPE_CHECKING, Literal, Union

import cftime
import pandas as pd
import xarray as xr

if TYPE_CHECKING:
    from wam2layers.config import Config

# Define a type alias for easier reuse
CfDateTime = Union[
    cftime.DatetimeNoLeap,
    cftime.Datetime360Day,
    cftime.DatetimeAllLeap,
    cftime.DatetimeGregorian,
    cftime.DatetimeJulian,
    cftime.DatetimeProlepticGregorian,
]


def subclass_cftime(dt: cftime.datetime) -> CfDateTime:
    """Convert the base cftime.datetime to the correct subclass.

    Base cftime.datetime objects don't place nice with xarray, by converting them to
    the specific subclasses things work better.

    This code is obsolete once https://github.com/pydata/xarray/pull/8942 is merged.
    """
    ordial = dt.toordinal(fractional=True)  # julian day ordinal
    if dt.calendar in ["noleap", "365_day"]:
        return cftime.DatetimeNoLeap.fromordinal(ordial, dt.calendar, dt.has_year_zero)
    elif dt.calendar == "360_day":
        return cftime.Datetime360Day.fromordinal(ordial, dt.calendar, dt.has_year_zero)
    elif dt.calendar in ["all_leap", "366_day"]:
        return cftime.DatetimeAllLeap.fromordinal(ordial, dt.calendar, dt.has_year_zero)
    elif dt.calendar in ["standard", "gregorian"]:
        return cftime.DatetimeGregorian.fromordinal(
            ordial, dt.calendar, dt.has_year_zero
        )
    elif dt.calendar == "julian":
        return cftime.DatetimeJulian.fromordinal(ordial, dt.calendar, dt.has_year_zero)
    elif dt.calendar == "proleptic_gregorian":
        return cftime.DatetimeProlepticGregorian.fromordinal(
            ordial, dt.calendar, dt.has_year_zero
        )
    else:
        raise ValueError


def has_cftime(ds: xr.Dataset) -> bool:
    """Check if the dataset has dimension 'time' with cftime values."""
    return "cftime" in str(type(ds["time"].values[0]))


def get_calendar_type(ds: xr.Dataset) -> str:
    """Look up the calendar type using cftime."""
    if has_cftime(ds):
        return ds["time"].values[0].calendar
    else:
        return "standard"


def validate_calendar_type(config: "Config") -> None:
    file = template_to_files(config, var="ua")[0]
    data_calendar = get_calendar_type(xr.open_dataset(file))
    if data_calendar != config.calendar:
        msg = (
            "Calendar mismatch detected:\n"
            f"    The configuration file has calendar type '{config.calendar}', "
            f"while the data has calendar '{data_calendar}'.\n"
            "    This can lead to unpredictable errors!"
        )
        warnings.warn(msg)


def template_to_files(config: "Config", var: str) -> list[Path]:
    template = Path(config.filename_template.format(variable=var))
    return list(template.parent.glob(template.name))


def fix_deprecated_frequency(freq: str) -> str:
    """Auto-correct deprecated pandas frequencies for backward compatibility.

    Pandas some time ago deprecated some frequencies:

    > Deprecated since version 2.2.0: Aliases H, BH, CBH, T, S, L, U, and N are
    > deprecated in favour of the aliases h, bh, cbh, min, s, ms, us, and ns.

    This fix recognizes the old format and converts it to the new one.
    """
    deprecated_map = {
        r"(?i)^(\d*)CBH$": r"\1cbh",
        r"(?i)^(\d*)BH$": r"\1bh",
        r"(?i)^(\d*)H$": r"\1h",
        r"(?i)^(\d*)T$": r"\1min",
        r"(?i)^(\d*)S$": r"\1s",
        r"(?i)^(\d*)L$": r"\1ms",
        r"(?i)^(\d*)U$": r"\1us",
        r"(?i)^(\d*)N$": r"\1ns",
    }

    for pattern, repl in deprecated_map.items():
        if re.fullmatch(pattern, freq):
            return re.sub(pattern, repl, freq)
    return freq


def round_cftime(
    date: CfDateTime, freq: str, how: Literal["nearest", "floor", "ceil"]
) -> CfDateTime:
    """Round a CF datetime index to a nearest frequency."""
    if how == "nearest":
        return xr.CFTimeIndex([date], calendar=date.calendar).round(freq)[0]
    if how == "ceil":
        return xr.CFTimeIndex([date], calendar=date.calendar).ceil(freq)[0]
    if how == "floor":
        return xr.CFTimeIndex([date], calendar=date.calendar).floor(freq)[0]


def cftime_from_iso(timestamp: str, calendar: str) -> CfDateTime:
    if timestamp.count(":") == 1:
        fmt = "%Y-%m-%dT%H:%M"
    else:
        fmt = "%Y-%m-%dT%H:%M:%S"
    return subclass_cftime(cftime.datetime.strptime(timestamp, fmt, calendar))


def cftime_from_timestamp(timestamp: pd.Timestamp, calendar: str) -> CfDateTime:
    return cftime_from_iso(timestamp.isoformat(), calendar)
