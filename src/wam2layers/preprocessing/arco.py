"""Load ARCO-ERA5: preprocessing from zarr."""

from time import time

import numpy as np
import pandas as pd
import xarray as xr

from wam2layers.config import Config
from wam2layers.preprocessing.cmip import LATENT_HEAT
from wam2layers.preprocessing.utils import calculate_humidity

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
    "surface_latent_heat_flux",
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


def get_input_data(datetime: pd.Timestamp, config: Config):
    """Get the preprocessing input data for ARCO-ERA5."""
    t0 = time()
    assert config.level_type == "pressure_levels"

    data = xr.open_zarr(
        "gs://gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3",
        chunks=None,
        storage_options=dict(token="anon"),
    )

    data = data.sel(time=datetime, method="nearest")

    ## Note: it is possible to reduce the domain, but due to the chunking (all lats/lons
    ## are in a single chunk, it won't save you any data transfer time). If storage
    ## allows: do the whole world.
    # if config.tracking_domain is not None:  # Reduce spatial domain
    #     data=data.sel(
    #         latitude=slice(config.tracking_domain.north, config.tracking_domain.south),
    #         longitude=slice(config.tracking_domain.west, config.tracking_domain.east),
    #     )
    data = data.sel(level=config.levels)  # Only select wanted levels

    data = data[REQUIRED_VARS]
    data["qs"] = calculate_humidity(
        data["2m_dewpoint_temperature"], data["surface_pressure"]
    )

    evap = data["surface_latent_heat_flux"] / LATENT_HEAT
    data["precip"] = np.maximum(data["total_precipitation"], 0) + np.maximum(evap, 0)
    data["evap"] = np.abs(np.minimum(evap, 0))

    data["level"] = data["level"] * 100  # hPa --> Pa
    data["level"].attrs.update(units="Pa")

    data = data.rename_vars(ARCO2WAM)
    data = data.drop_vars(
        [
            "2m_temperature",
            "2m_dewpoint_temperature",
            "surface_latent_heat_flux",
            "total_precipitation",
        ]
    )
    data.load()  # Load ARCO data into memory
    print(f"Elapsed time ARCO load: {time() - t0:.1f} s")
    return data, {}
