import xarray as xr


def load_region(config):
    # TODO: make variable name more generic
    return xr.open_dataset(config.region).source_region
