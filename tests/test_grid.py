import numpy as np
import xarray as xr

from wam2layers.utils.grid import get_boundary


def test_get_boundary():
    a = np.arange(12).reshape(3, 4)
    a_xr = xr.DataArray(a)
    result = get_boundary(a_xr)
    expected = np.array(
        [
            [True, True, True, True],
            [True, False, False, True],
            [True, True, True, True],
        ]
    )
    np.testing.assert_array_equal(result, expected)


def test_get_boundary_periodic():
    a = np.arange(12).reshape(3, 4)
    a_xr = xr.DataArray(a)
    result = get_boundary(a_xr, periodic=True)
    expected = np.array(
        [
            [True, True, True, True],
            [False, False, False, False],
            [True, True, True, True],
        ]
    )
    np.testing.assert_array_equal(result, expected)
