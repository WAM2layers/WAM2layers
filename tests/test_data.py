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
pressure_4d = dummy_data((5, 4, 3, 2), monotonic_decreasing=True)
temperature_4d = dummy_data((5, 4, 3, 2))
assert np.all(pressure_4d.diff("lev") <= 0), "Pressure levels should be monotonic"
