"""Load ARCO-ERA5: preprocessing from zarr."""

import datetime as pydt
from functools import cache

import numpy as np
import xarray as xr

from wam2layers import __version__
from wam2layers.config import Config
from wam2layers.preprocessing.utils import accumulation_to_flux, calculate_humidity
from wam2layers.utils.calendar import CfDateTime

PREPROCESS_ATTRS = {
    "title": "ARCO-ERA5 data preprocessed for use in WAM2layers",
    "history": (
        f"created on {pydt.datetime.now(pydt.timezone.utc):%Y-%m-%dT%H:%M:%SZ} "
        f"using wam2layers version {__version__}."
    ),
    "source": (
        "ARCO-ERA5: An Analysis-Ready Cloud-Optimized Reanalysis Dataset, "
        "https://github.com/google-research/arco-era5, "
        "ECMWF Reanalysis v5 (ERA5), "
        "www.ecmwf.int/en/forecasts/dataset/ecmwf-reanalysis-v5"
    ),
    "references": (
        "doi.org/10.5281/zenodo.7010594, "
        "doi.org/10.5194/esd-5-471-2014, "
        "https://ams.confex.com/ams/103ANNUAL/meetingapp.cgi/Paper/415842"
    ),
    "Conventions": "CF-1.6",
}


REQUIRED_VARS = [
    "u_component_of_wind",
    "v_component_of_wind",
    "specific_humidity",
    "10m_u_component_of_wind",
    "10m_v_component_of_wind",
    "surface_pressure",
    "2m_temperature",
    "2m_dewpoint_temperature",
    "total_column_water",
    "total_precipitation",
    "evaporation",
]

ARCO2WAM = {
    "u_component_of_wind": "u",
    "v_component_of_wind": "v",
    "specific_humidity": "q",
    "10m_u_component_of_wind": "us",
    "10m_v_component_of_wind": "vs",
    "surface_pressure": "ps",
    "total_column_water": "tcw",
}


def get_input_data(datetime: CfDateTime, config: Config):
    """Get the preprocessing input data for ARCO-ERA5."""
    if config.level_type != "pressure_levels":
        msg = "Currently the ARCO-ERA5 preprocessor only supports pressure level data."
        raise ValueError(msg)

    if config.calendar != "proleptic_gregorian":
        msg = (
            "ARCO-ERA5 data is on a 'proleptic_gregorian' calendar.\n"
            "Please set this calendar in your configuration file."
        )
        raise ValueError(msg)

    ds = open_arco_pl()
    data = ds.sel(time=datetime).compute()

    data = data.sel(level=config.levels)  # Only select wanted levels

    data["qs"] = calculate_humidity(
        data["2m_dewpoint_temperature"], data["surface_pressure"]
    )

    # Transfer negative (originally positive) values of evap to precip
    data["precip"] = np.maximum(data["total_precipitation"], 0) + np.maximum(
        data["evaporation"], 0
    )

    # Change sign convention to all positive,
    data["evap"] = np.abs(np.minimum(data["evaporation"], 0))

    data["precip"] = accumulation_to_flux(
        data["precip"], input_frequency=config.input_frequency
    )
    data["evap"] = accumulation_to_flux(
        data["evap"], input_frequency=config.input_frequency
    )

    data["level"] = data["level"] * 100  # hPa --> Pa
    data["level"].attrs.update(units="Pa")

    data = data.rename_vars(ARCO2WAM)
    data = data.drop_vars(
        [
            "2m_temperature",
            "2m_dewpoint_temperature",
            "evaporation",
            "total_precipitation",
        ]
    )
    data.load()  # Load ARCO data into memory
    return data, PREPROCESS_ATTRS


@cache
def open_arco_pl():
    """Open the required variables from ARCO in an efficient manner."""
    data = xr.open_zarr(
        "gs://gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3",
        chunks=None,
        use_cftime=True,
        storage_options={"token": "anon"},
    )[REQUIRED_VARS]

    return data.sel(
        time=slice(data.attrs["valid_time_start"], data.attrs["valid_time_stop"])
    )
