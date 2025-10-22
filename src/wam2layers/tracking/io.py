import logging
import time
import traceback
from copy import deepcopy
from datetime import datetime as pydt
from functools import lru_cache
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import xarray as xr

from wam2layers import __version__
from wam2layers.config import Config
from wam2layers.preprocessing.utils import add_bounds
from wam2layers.reference.variables import DATA_ATTRIBUTES
from wam2layers.utils.calendar import CfDateTime, round_cftime
from wam2layers.utils.shapefile_mask import create_mask

logger = logging.getLogger(__name__)


def input_path(date: CfDateTime, input_dir: Path):
    return f"{input_dir}/{date.strftime('%Y-%m-%d')}_fluxes_storages.nc"


def output_path(date: CfDateTime, config: Config):
    output_dir = config.output_folder
    mode = "backtrack" if config.tracking_direction == "backward" else "forwardtrack"
    return f"{output_dir}/{mode}_{date.strftime('%Y-%m-%dT%H-%M')}.nc"


def fully_load_and_detach(ds: xr.Dataset) -> xr.Dataset:
    """Force all variables into memory and detach from OPeNDAP/netCDF backend."""
    # Extract data variables with explicit dimension names
    loaded_vars = {var: (ds[var].dims, ds[var].values) for var in ds.data_vars}

    # Extract coords similarly (some coords are multi-dimensional)
    loaded_coords = {
        coord: (ds.coords[coord].dims, ds.coords[coord].values) for coord in ds.coords
    }
    return xr.Dataset(data_vars=loaded_vars, coords=loaded_coords, attrs=ds.attrs)


def select_subdomain(ds: xr.Dataset, bbox: str):
    """Limit the data to a subdomain."""
    west, south, east, north = [float(p) for p in bbox.split(", ")]
    if west < east:
        return ds.sel(
            longitude=slice(west, east),  # -180 to 180
            latitude=slice(north, south),  # 90 - -90 (decreasing)
            drop=True,
        )
    else:
        # Need to roll coordinates to maintain a continuous longitude
        # Roll the first lon that's NOT in the bbox to the front of the array
        # That way, all values that ARE in the bbox, must be contiguous.
        in_bbox = (ds.longitude < east) | (ds.longitude > west)
        n_roll = in_bbox.argmin().item()  # first false value
        ds_rolled = ds.roll(longitude=-n_roll, roll_coords=True)
        in_bbox = (ds_rolled.longitude < east) | (ds_rolled.longitude > west)
        return ds_rolled.sel(longitude=in_bbox, latitude=slice(north, south))


def read_or_retry(filename: str) -> xr.Dataset:
    """Read a file, retrying if it fails.

    Avoids spurious errors in decoding time when using OpenDAP

    Retry behaviour copied from xarray.backends.common.robust_getitem
    """
    max_retries = 8  # ~ 4 minutes maximum, so far always worked within 6 tries.
    initial_delay = 500  # milliseconds

    for n in range(max_retries + 1):
        try:
            return xr.open_dataset(
                filename, engine="netcdf4", cache=False, use_cftime=True
            )
        except IndexError:
            if n == max_retries:
                raise
            base_delay = initial_delay * 2**n
            next_delay = base_delay + np.random.randint(base_delay)
            msg = (
                f"file read failed, waiting {next_delay} ms before trying again "
                f"({max_retries - n} tries remaining). {traceback.format_exc(limit=2, chain=False)}"
            )
            logger.debug(msg)
            time.sleep(1e-3 * next_delay)


@lru_cache(maxsize=6)
def _load_slice_with_cache(filename, time, subset: str, bbox: Optional[str] = None):
    """Load data slice from file."""
    variables = {
        "fluxes": ["fx_upper", "fx_lower", "fy_upper", "fy_lower", "evap", "precip"],
        "states": ["s_upper", "s_lower"],
    }[subset]

    logger.debug(f"Loading {subset} at time {time} from {filename}")

    ds = read_or_retry(filename)
    ds = ds.sel(time=time).drop("time").squeeze()
    ds = ds[variables]

    if bbox is not None:
        ds = select_subdomain(ds, bbox)

    return fully_load_and_detach(ds)


def load_slice(t: CfDateTime, subset: str, config: Config) -> xr.Dataset:
    filename = input_path(t, config.preprocessed_data_folder)
    time = t.strftime("%Y-%m-%dT%H:%M:%S")
    domain = str(config.tracking_domain) if config.tracking_domain else None

    # Avoid modifying the cached slice!
    slice = _load_slice_with_cache(filename, time, subset, domain)
    return deepcopy(slice)


def load_data(t: CfDateTime, config: Config, subset: str = "fluxes") -> xr.Dataset:
    """Load variable at t, interpolate if needed."""
    t0 = round_cftime(t, config.input_frequency, "floor")
    t1 = round_cftime(t, config.input_frequency, "ceil")

    if t == t0 or t == t1:
        return load_slice(t, subset, config)

    da0 = load_slice(t0, subset, config)
    da1 = load_slice(t1, subset, config)
    return da0 + (t - t0) / (t1 - t0) * (da1 - da0)


def load_tagging_region(config: Config, t: Optional[CfDateTime] = None):
    if config.tagging_region.suffix == ".shp":
        logger.info("Generating mask from shapefile...")
        ds = xr.open_dataset(next(config.preprocessed_data_folder.glob("*.nc")))
        tagging_region = create_mask(ds, config.tagging_region)
        logger.info("Writing mask to debug folder.")
        (config.output_folder / "debug").mkdir(exist_ok=True)
        tagging_region.to_dataset().to_netcdf(
            config.output_folder / "debug" / "tagging_region.nc"
        )
    else:
        tagging_region = xr.open_dataarray(config.tagging_region)

    if config.tracking_domain is not None:
        tagging_region = select_subdomain(tagging_region, str(config.tracking_domain))

    if t is not None:
        return tagging_region.sel(time=t, method="nearest")
    return tagging_region


def add_time_bounds(ds: xr.Dataset, t: CfDateTime, config: Config):
    """Compute the time bounds and add them to the output dataset."""
    dt = pd.Timedelta(config.output_frequency)
    if config.tracking_direction == "forward":
        ds["time_bnds"] = (("time", "bnds"), np.array(((t - dt, t),)))
    else:
        ds["time_bnds"] = (("time", "bnds"), np.array(((t + dt, t),)))


def get_main_attrs(attrs: dict, config: Config):
    """Create the main dataset attributes, including appending to the history."""
    new_history = (
        f"created on {pydt.utcnow():%Y-%m-%dT%H:%M:%SZ} "
        f"using wam2layers version {__version__}; "
    )
    new_attrs = {
        "title": f"WAM2Layers {config.tracking_direction}tracking output file",
        "history": (
            new_history + attrs["history"] if "history" in attrs else new_history
        ),
        "source": (
            "Moisture tracking applied to original data; "
            f"{attrs['source'] if 'source' in attrs else 'UNKNOWN'}"
        ),
        "references": "doi.org/10.5281/zenodo.7010594, doi.org/10.5194/esd-5-471-2014",
        "Conventions": "CF-1.6",
        "WAM2Layers_config": config.to_string(),
    }

    return new_attrs


def write_output(
    output: xr.Dataset,
    t: CfDateTime,
    config: Config,
) -> None:
    # TODO: add back (and cleanup) coordinates and units
    path = output_path(t, config)
    logger.info(f"{t} - Writing output to file {path}")

    output = output.assign_coords({"time": t})
    output = output.expand_dims(dim="time")

    for var in output.data_vars:
        output[var].attrs.update(DATA_ATTRIBUTES[str(var)])

    add_bounds(output)
    add_time_bounds(output, t, config)
    output.attrs = get_main_attrs(output.attrs, config)

    comp = dict(zlib=True, complevel=8)
    encoding = {var: comp for var in output.data_vars}
    time_encoding = {"units": "seconds since 1900-01-01"}  # TODO: check cftime output
    encoding["time"] = time_encoding
    encoding["time_bnds"] = time_encoding

    for var in output.data_vars:
        if var != "time_bnds":
            output[var] = output[var].astype("float32")
    output.to_netcdf(path, encoding=encoding)
