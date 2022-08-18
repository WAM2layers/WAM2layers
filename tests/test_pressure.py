import xarray as xr
import numpy as np


def dummy_data(shape, fill_value=None, monotonic_decreasing=False):
    """Create dummy data with a given shape and fill value."""
    if fill_value is not None:
        data = np.zeros(shape) + fill_value
    else:
        data = np.random.randint(0, 9, shape)

    if monotonic_decreasing:
        data = np.flip(np.sort(data, axis=1), axis=1)

    if len(shape) == 4:
        return xr.DataArray(data, dims=["time", "lev", "lat", "lon"], coords=[range(i) for i in shape])
    elif len(shape) == 3:
        return xr.DataArray(data, dims=["time", "lat", "lon"], coords=[range(i) for i in shape])

# Make some dummy data to play with

# 4d arrays valid at pressure levels
pressure_4d = dummy_data((5, 4, 3, 2), monotonic_decreasing=True)
temperature_4d = dummy_data((5, 4, 3, 2))
assert np.all(pressure_4d.diff("lev") <= 0), "Pressure levels should be monotonic"

# 3d Arrays valid at the surface
pressure_surface = dummy_data((5, 3, 2), fill_value=11)
temperature_surface = dummy_data((5, 3, 2), fill_value=-5)
pressure_surface["lev"] = -1  # note: lev is just counting levels, no physical meaning
temperature_surface["lev"] = -1

# 3d arrays valid at top of atmosphere
pressure_toa = dummy_data((5, 3, 2), fill_value=-1)
temperature_toa = temperature_4d.isel(lev=-1)
pressure_toa["lev"] = 5  # note: lev is just counting levels, no physical meaning
temperature_toa["lev"] = 5

# Create new 4d array combining surface, pressure levels and top of atmosphere
# Beware: this doesn't maintain dimension order, level is the last dimension now!
pressure_4d = xr.concat([pressure_surface, pressure_4d, pressure_toa], dim="lev")
temperature_4d = xr.concat([temperature_surface, temperature_4d, temperature_toa], dim="lev")

# Sort the data
# TODO: with some tricks this can be done in xarray directly
sorted_p, sorted_t = sortby_ndarray(pressure_4d.values, temperature_4d.values, axis=3)
pressure_4d.values = sorted_p[..., ::-1]  # restore descending order
temperature_4d.values = sorted_p[..., ::-1]

pressure_boundary = dummy_data((5, 3, 2), fill_value=5)
temperature_boundary = dummy_data((5, 3, 2), fill_value=0)  # fill with dummy data for now, then interpolate to set the values
pressure_boundary["lev"] = -2  # note: will insert at bottom first, sort later
temperature_boundary["lev"] = -2  # note: will insert at bottom first, sort later
temperature_boundary.values = interpolate(pressure_boundary.values, pressure_4d.values, temperature_4d.values, axis=3, descending=True)

# Create new 4d array adding the boundary to the already existing levels
# note: will insert at bottom first, sort later
pressure_4d = xr.concat([pressure_boundary, pressure_4d], dim="lev")
temperature_4d = xr.concat([temperature_boundary, temperature_4d], dim="lev")

# Sort the data (again..)
# TODO: with some tricks this can be done in xarray directly
sorted_p, sorted_t = sortby_ndarray(pressure_4d.values, temperature_4d.values, axis=3)
pressure_4d.values = sorted_p[..., ::-1]  # restore descending order
temperature_4d.values = sorted_p[..., ::-1]

assert np.all(pressure_4d.diff("lev") <= 0), "Pressure levels should be monotonic"
