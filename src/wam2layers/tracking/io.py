import logging

import xarray as xr

logger = logging.getLogger(__name__)


def input_path(date, config):
    input_dir = config.preprocessed_data_folder
    return f"{input_dir}/{date.strftime('%Y-%m-%d')}_fluxes_storages.nc"


def output_path(date, config):
    output_dir = config.output_folder
    return f"{output_dir}/backtrack_{date.strftime('%Y-%m-%dT%H-%M')}.nc"


def read_data_at_date(d, config):
    """Load input data for given date."""
    file = input_path(d, config)
    return xr.open_dataset(file, cache=False)


def read_data_at_time(t, config):
    """Get a single time slice from input data at time t."""
    ds = read_data_at_date(t, config)
    return ds.sel(time=t, drop=True)


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


def load_region(config):
    # TODO: make variable name more generic
    return xr.open_dataset(config.region).tagging_region


def write_output(output, t, config):
    # TODO: add back (and cleanup) coordinates and units
    path = output_path(t, config)
    logger.info(f"{t} - Writing output to file {path}")
    output.astype("float32").to_netcdf(path)

    # Flush previous outputs
    output[["e_track", "tagged_precip"]] *= 0
