from functools import cache
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from wam2layers.config import Config
from wam2layers.utils.calendar import CfDateTime, template_to_files

WAM2CMIP = {  # wam2layers name -> cmip name
    "q": "hus",
    "e": "hfls",
    "tp": "pr",
    "ps": "ps",
    "u": "ua",
    "v": "va",
    "qs": "huss",
    "us": "uas",
    "vs": "vas",
}
CMIP2WAM = {var: key for key, var in WAM2CMIP.items()}  # cmip name --> wam2layers name

COORD_NAMES = {
    "lat": "latitude",
    "lon": "longitude",
    "plev": "level",
}

THREE_HR_VARS = ("pr", "hfls")

LATENT_HEAT = 2.45e6  # J/kg (at ~20 degC, ranges from 2.5 (0 degC) to 2.4 (50 degC)
# 1.91846e6*((temp/(temp-33.91))**2)


def get_input_data(datetime: CfDateTime, config: Config):
    assert config.level_type == "pressure_levels"

    data = xr.Dataset()
    q = load_cmip_var("q", datetime, config)  # in kg kg-1
    u = load_cmip_var("u", datetime, config)  # in m/s
    v = load_cmip_var("v", datetime, config)  # in m/s
    levels = common_levels(q["level"], u["level"], v["level"])

    data["q"] = load_cmip_var("q", datetime, config).sel(level=levels)
    data["u"] = load_cmip_var("u", datetime, config).sel(level=levels)
    data["v"] = load_cmip_var("v", datetime, config).sel(level=levels)
    data["ps"] = load_cmip_var("ps", datetime, config)  # in Pa

    # Pressure level input vars
    data["qs"] = load_cmip_var("qs", datetime, config)  # kg kg-1
    data["us"] = load_cmip_var("us", datetime, config)  # in m/s
    data["vs"] = load_cmip_var("vs", datetime, config)  # in m/s

    data["precip"] = load_cmip_var("tp", datetime, config)
    data["evap"] = load_cmip_var("e", datetime, config)

    data["precip"] = np.maximum(data["precip"], 0) + np.minimum(data["evap"], 0)
    data["evap"] = np.abs(np.maximum(data["evap"], 0))

    return data


@cache  # Cache this function: keep output in memory for reuse
def open_var(files: tuple[Path], var: str) -> xr.DataArray:
    """Open a certain CMIP variable, which can be split over multiple files time-wise.

    To speed up this function, the output is cached so it can be reused
    in the following timesteps.

    Args:
        files: The netCDF files corresponding to this variable.
        var: The name of the varialbe in the netCDF file.
    """
    datasets = [
        xr.open_dataset(
            file,
            chunks={"time": 1, "latitude": -1, "longitude": -1},
            use_cftime=True,
        )
        for file in files
    ]

    ds = xr.concat(datasets, dim="time") if len(datasets) > 0 else datasets[0]
    return ds[var]


def load_cmip_var(variable: str, datetime: CfDateTime, config: Config):
    """Load a single CMIP input data variable.

    Args:
        variable: Variable (wam2layers) name
        datetime: Time at which you want the variable
        config: WAM2layers configuration

    Returns:
        DataArray of the variable at the specified time
    """
    var = WAM2CMIP[variable]
    input_files = tuple(template_to_files(config, var))

    data = open_var(input_files, var)

    if var in THREE_HR_VARS:
        data = resample_3hr_var(data)
    if var == "hfls":
        data = data / LATENT_HEAT

    data = data.rename(CMIP2WAM[data.name])
    for old_name, new_name in COORD_NAMES.items():
        if old_name in data.coords:
            data = data.rename({old_name: new_name})

    return data.sel(time=datetime)


def resample_3hr_var(data: xr.DataArray) -> xr.DataArray:
    """Resample 3-hour flux variables to 6 hour intervals."""
    assert_correct_starttime(data)
    data = data.copy(deep=True)

    data["time"] = np.repeat(
        data["time"].isel(time=slice(None, None, 2)) + pd.Timedelta("4h30m"), repeats=2
    )

    return (data.isel(time=slice(0, None, 2)) + data.isel(time=slice(1, None, 2))) / 2


def assert_correct_starttime(data: xr.DataArray):
    if (data["time"][0].dt.hour, data["time"][0].dt.minute) != (1, 30):
        msg = (
            f"Incorrect start time for 3 hour variable '{data.name}'.\n"
            f"First time value should start at 01:30:00, but is {data['time'][0]}"
        )
        raise ValueError(msg)


def common_levels(*levels: xr.DataArray) -> np.ndarray:
    """Find common pressure/model levels between multiple variables."""
    lvs = iter(levels)
    common = set(next(lvs).to_numpy())
    for lv in lvs:
        common = common.intersection(lv.to_numpy())
    return np.array(list(common))
