"""Loading ERA5 data for preprocessing."""
from datetime import datetime as pydt
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from wam2layers import __version__
from wam2layers.config import Config
from wam2layers.preprocessing import logger
from wam2layers.preprocessing.utils import (
    accumulation_to_flux,
    calculate_humidity,
    log_once,
    midpoints,
)

PREPROCESS_ATTRS = {
    "title": "ERA5 data preprocessed for use in WAM2layers",
    "history": (
        f"created on {pydt.utcnow():%Y-%m-%dT%H:%M:%SZ} "
        f"using wam2layers version {__version__}."
    ),
    "source": (
        "ECMWF Reanalysis v5 (ERA5), "
        "www.ecmwf.int/en/forecasts/dataset/ecmwf-reanalysis-v5"
    ),
    "references": "doi.org/10.5281/zenodo.7010594, doi.org/10.5194/esd-5-471-2014",
    "Conventions": "CF-1.6",
}


def get_input_data(datetime: pd.Timestamp, config: Config):
    """Get the preprocessing input data for ERA5."""
    data = xr.Dataset()
    data["q"] = load_data("q", datetime, config)  # in kg kg-1
    data["u"] = load_data("u", datetime, config)  # in m/s
    data["v"] = load_data("v", datetime, config)  # in m/s
    data["ps"] = load_data("sp", datetime, config)  # in Pa
    if config.level_type == "pressure_levels":
        d2m = load_data("d2m", datetime, config)  # Dew point in K
        data["qs"] = calculate_humidity(d2m, data["ps"])  # kg kg-1
        data["us"] = load_data("u10", datetime, config)  # in m/s
        data["vs"] = load_data("v10", datetime, config)  # in m/s
    if config.level_type == "model_levels":
        data["dp"] = get_dp_modellevels(data["ps"], config.levels)

    data["precip"], data["evap"] = preprocess_precip_and_evap(datetime, config)
    try:
        data["tcw"] = load_data("tcw", datetime, config)
    except FileNotFoundError:
        log_once(
            logger,
            "Total column water input data not available; using water vapour only",
        )

    if config.level_type == "pressure_levels":
        data["level"] = data["level"] * 100  # hPa --> Pa
        data["level"].attrs.update(units="Pa")

    return data, PREPROCESS_ATTRS


def load_data(variable: str, datetime: pd.Timestamp, config: Config) -> xr.DataArray:
    """Load data for given variable and date."""
    template = config.filename_template

    # If it's a 4d variable, we need to set the level type
    if variable in ["u", "v", "q"]:
        if config.level_type == "model_levels":
            prefix = "_ml"
        elif config.level_type == "pressure_levels":
            prefix = "_pl"
    # If it's a 3d variable we do not need to set the level type
    else:
        prefix = ""

    # Load data
    filepath = template.format(
        year=datetime.year,
        month=datetime.month,
        day=datetime.day,
        hour=datetime.hour,
        minute=datetime.minute,
        levtype=prefix,
        variable=variable,
    )
    da = xr.open_dataarray(filepath).sel(time=datetime.strftime("%Y%m%d%H%M"))

    if "lev" in da.coords:
        da = da.rename(lev="level")

    # If it's 4d data we want to select a subset of the levels
    if variable in ["u", "v", "q"] and isinstance(config.levels, list):
        return da.sel(level=config.levels)

    return da


def preprocess_precip_and_evap(
    datetime: pd.Timestamp, config: Config
) -> tuple[xr.DataArray, xr.DataArray]:
    """Load and pre-process precipitation and evaporation."""
    # All incoming units are accumulations (in m) since previous time step
    evap = load_data("e", datetime, config)
    precip = load_data("tp", datetime, config)  # total precipitation

    # Transfer negative (originally positive) values of evap to precip
    precip = np.maximum(precip, 0) + np.maximum(evap, 0)
    precip = precip

    # Change sign convention to all positive,
    evap = np.abs(np.minimum(evap, 0))

    precip = accumulation_to_flux(precip, input_frequency=config.input_frequency)
    evap = accumulation_to_flux(evap, input_frequency=config.input_frequency)
    return precip, evap


def load_era5_ab():
    """Load a and b coefficients."""
    table = Path(__file__).parent / "tableERA5model_to_pressure.csv"
    df = pd.read_csv(table)
    a = df["a [Pa]"].values
    b = df.b.values
    return a, b


def get_edges(era5_modellevels):
    """Get the values of a and b at the edges of a subset of ERA5 modellevels."""
    a, b = load_era5_ab()

    if era5_modellevels == "all":
        return a, b

    # With python indexing starting at 0 and 137 full levels,
    # count from 0-136 instead of 1-137.
    era5_modellevels = [i - 1 for i in era5_modellevels]

    # Calculate the values of a and be at midpoints and extract levels
    a_full = midpoints(a)[era5_modellevels]
    b_full = midpoints(b)[era5_modellevels]

    # Calculate the values of a and b at the edges *between* these levels
    a_edge = midpoints(a_full)
    b_edge = midpoints(b_full)

    # Reinsert the original outer edges
    a_edge = np.block([a[0], a_edge, a[-1]])
    b_edge = np.block([b[0], b_edge, b[-1]])

    return a_edge, b_edge


def get_dp_modellevels(sp, levels):
    """Calculate pressure jump over subset of era5 modellevels."""
    a, b = get_edges(levels)

    # Use dummy level -1 to avoid collision with existing levels
    # first coord will be discarded anyway on taking the diff below.
    a = xr.DataArray(a, coords={"level": np.insert(levels, 0, -1)})
    b = xr.DataArray(b, coords={"level": np.insert(levels, 0, -1)})

    p_edge = a + b * sp

    # Calculate the difference between the pressure levels
    dp = p_edge.diff("level")  # in Pa

    return dp
