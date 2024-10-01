import numpy as np
import pandas as pd
import xarray as xr

from wam2layers.config import Config

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

SECONDS_6HR = 6 * 3600
LATENT_HEAT = 2.45e6  # J/kg (at ~20 degC, ranges from 2.5 (0 degC) to 2.4 (50 degC)
# 1.91846e6*((temp/(temp-33.91))**2)


def get_input_data(datetime, config):
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

    data["precip"] = np.maximum(data["precip"], 0) + np.maximum(data["evap"], 0)
    data["evap"] = np.abs(np.minimum(data["evap"], 0))

    return data


def load_cmip_var(variable, datetime: pd.Timestamp, config: Config):
    """Load a single CMIP input data variable.

    Args:
        variable: Variable (wam2layers) name
        datetime: Time at which you want the variable
        config: WAM2layers configuration

    Returns:
        DataArray of the variable at the specified time
    """
    var = WAM2CMIP[variable]
    fname_pattern = config.filename_template.format(variable=var)

    ds = xr.open_mfdataset(fname_pattern, chunks="auto")
    data = ds[var]

    if var in THREE_HR_VARS:
        data = accumulate_3hr_var(data)
    if var == "hfls":
        data = data / LATENT_HEAT

    data = data.rename(CMIP2WAM[data.name])
    for old_name, new_name in COORD_NAMES.items():
        if old_name in data.coords:
            data = data.rename({old_name: new_name})

    return data.sel(time=datetime.strftime("%Y%m%d%H%M"))


def accumulate_3hr_var(data: xr.DataArray) -> xr.DataArray:
    """Accumulates 3-hour variables to 6 hour intervals."""
    assert_correct_starttime(data)
    data = data.copy(deep=True)

    data["time"] = np.repeat(
        data["time"].isel(time=slice(None, None, 2)) + pd.Timedelta("4h30m"), repeats=2
    )

    return (
        (data.isel(time=slice(0, None, 2)) + data.isel(time=slice(1, None, 2)))
        / 2
        * SECONDS_6HR
    )


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
