import logging
from datetime import datetime as pydt

import numpy as np
import pandas as pd
import xarray as xr

from wam2layers import __version__
from wam2layers.config import Config
from wam2layers.preprocessing.utils import add_bounds
from wam2layers.reference.variables import DATA_ATTRIBUTES

logger = logging.getLogger(__name__)


def input_path(date, config):
    input_dir = config.preprocessed_data_folder
    return f"{input_dir}/{date.strftime('%Y-%m-%d')}_fluxes_storages.nc"


def output_path(date, config):
    output_dir = config.output_folder
    mode = "backtrack" if config.tracking_direction == "backward" else "forwardtrack"
    return f"{output_dir}/{mode}_{date.strftime('%Y-%m-%dT%H-%M')}.nc"


def read_data_at_date(d, config):
    """Load input data for given date."""
    file = input_path(d, config)
    ds = xr.open_dataset(file, cache=False)
    if config.tracking_domain is not None:
        return select_subdomain(ds, config.tracking_domain)
    return ds


def read_data_at_time(t, config):
    """Get a single time slice from input data at time t."""
    ds = read_data_at_date(t, config)
    return ds.sel(time=t, drop=True)


def select_subdomain(ds, bbox):
    """Limit the data to a subdomain."""
    if bbox.west < bbox.east:
        return ds.sel(
            longitude=slice(bbox.west, bbox.east),  # -180 to 180
            latitude=slice(bbox.north, bbox.south),  # 90 - -90 (decreasing)
            drop=True,
        )
    else:
        # Need to roll coordinates to maintain a continuous longitude
        # Roll the first lon that's NOT in the bbox to the front of the array
        # That way, all values that ARE in the bbox, must be contiguous.
        in_bbox = (ds.longitude < bbox.east) | (ds.longitude > bbox.west)
        n_roll = in_bbox.argmin().item()  # first false value
        ds_rolled = ds.roll(longitude=-n_roll, roll_coords=True)
        in_bbox = (ds_rolled.longitude < bbox.east) | (ds_rolled.longitude > bbox.west)
        return ds_rolled.sel(longitude=in_bbox, latitude=slice(bbox.north, bbox.south))


def load_data(t, config, subset="fluxes"):
    """Load variable at t, interpolate if needed."""
    variables = {
        "fluxes": ["fx_upper", "fx_lower", "fy_upper", "fy_lower", "evap", "precip"],
        "states": ["s_upper", "s_lower"],
    }

    t1 = t.ceil(config.input_frequency)
    da1 = read_data_at_time(t1, config)[variables[subset]]
    if t == t1:
        # this saves a lot of work if times are aligned with input
        return da1

    t0 = t.floor(config.input_frequency)
    da0 = read_data_at_time(t0, config)[variables[subset]]
    if t == t0:
        return da0

    return da0 + (t - t0) / (t1 - t0) * (da1 - da0)


def load_tagging_region(config, t=None):
    tagging_region = xr.open_dataarray(config.tagging_region)

    if config.tracking_domain is not None:
        tagging_region = select_subdomain(tagging_region, config.tracking_domain)

    if t is not None:
        return tagging_region.sel(time=t, method="nearest")
    return tagging_region


def add_time_bounds(ds: xr.Dataset, t: pd.Timestamp, config: Config):
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
    return {
        "title": f"WAM2Layers {config.tracking_direction}tracking output file",
        "history": new_history + attrs["history"]
        if "history" in attrs
        else new_history,
        "source": (
            "Moisture tracking applied to ERA5 dataset: ECMWF Reanalysis v5 (ERA5), "
            "www.ecmwf.int/en/forecasts/dataset/ecmwf-reanalysis-v5"
        ),
        "references": "doi.org/10.5281/zenodo.7010594, doi.org/10.5194/esd-5-471-2014",
        "Conventions": "CF-1.6",
        "WAM2Layers_config": config.model_dump_json(),
    }


def write_output(
    output: xr.Dataset,
    t: pd.Timestamp,
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
    time_encoding = {"units": "seconds since 1900-01-01"}
    encoding["time"] = time_encoding
    encoding["time_bnds"] = time_encoding

    for var in output.data_vars:
        if var != "time_bnds":
            output[var] = output[var].astype("float32")
    output.to_netcdf(path, encoding=encoding)
