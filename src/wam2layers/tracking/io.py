from functools import lru_cache

import xarray as xr

from wam2layers.tracking.backtrack import config


def input_path(date, config):
    input_dir = config.preprocessed_data_folder
    return f"{input_dir}/{date.strftime('%Y-%m-%d')}_fluxes_storages.nc"


def output_path(date, config):
    output_dir = config.output_folder
    return f"{output_dir}/backtrack_{date.strftime('%Y-%m-%dT%H-%M')}.nc"


# LRU Cache keeps the file open so we save a bit on I/O
@lru_cache(maxsize=2)
def read_data_at_date(d):
    """Load input data for given date."""
    file = input_path(d, config)
    return xr.open_dataset(file, cache=False)


# This one can't be cached as we'll be overwriting the content.
def read_data_at_time(t):
    """Get a single time slice from input data at time t."""
    ds = read_data_at_date(t)
    return ds.sel(time=t, drop=True)


def load_data(t, subset="fluxes"):
    """Load variable at t, interpolate if needed."""
    variables = {
        "fluxes": ["fx_upper", "fx_lower", "fy_upper", "fy_lower", "evap", "precip"],
        "states": ["s_upper", "s_lower"],
    }

    t1 = t.ceil(config.input_frequency)
    da1 = read_data_at_time(t1)[variables[subset]]
    if t == t1:
        # this saves a lot of work if times are aligned with input
        return da1

    t0 = t.floor(config.input_frequency)
    da0 = read_data_at_time(t0)[variables[subset]]
    if t == t0:
        return da0

    return da0 + (t - t0) / (t1 - t0) * (da1 - da0)


def load_region(config):
    # TODO: make variable name more generic
    return xr.open_dataset(config.region).source_region
