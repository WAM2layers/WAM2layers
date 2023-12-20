"""
Utility for appending to netcdf file along unlimited dimension.
Adapted from https://github.com/pydata/xarray/issues/1672#issuecomment-685222909
"""

import netCDF4
import xarray as xr


def _expand_variable(nc_variable, data, expanding_dim, nc_shape, added_size):
    # For time deltas, we must ensure that we use the same encoding as
    # what was previously stored.
    # We likely need to do this as well for variables that had custom
    # econdings too
    if hasattr(nc_variable, "calendar"):
        data.encoding = {
            "units": nc_variable.units,
            "calendar": nc_variable.calendar,
        }
    data_encoded = xr.conventions.encode_cf_variable(data)
    left_slices = data.dims.index(expanding_dim)
    right_slices = data.ndim - left_slices - 1
    nc_slice = (
        (slice(None),) * left_slices
        + (slice(nc_shape, nc_shape + added_size),)
        + (slice(None),) * (right_slices)
    )
    nc_variable[nc_slice] = data_encoded.data


def append_to_netcdf(filename, ds_to_append, expanding_dim="time"):
    with netCDF4.Dataset(filename, mode="a") as nc:
        nc_coord = nc[expanding_dim]
        nc_shape = len(nc_coord)
        added_size = len(ds_to_append[expanding_dim])
        variables, attrs = xr.conventions.encode_dataset_coordinates(ds_to_append)

        for name, data in variables.items():
            if expanding_dim not in data.dims:
                # Nothing to do, data assumed to the identical
                continue

            nc_variable = nc[name]
            _expand_variable(nc_variable, data, expanding_dim, nc_shape, added_size)
