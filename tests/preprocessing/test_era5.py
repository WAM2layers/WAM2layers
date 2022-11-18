from wam2layers.preprocessing.era5 import midpoints, get_edges, load_era5_ab
from numpy.testing import assert_array_equal, assert_almost_equal
import numpy as np
import pytest

def test_midponits():
    a = np.array([1, 2, 3])
    result = midpoints(a)
    expected = [1.5, 2.5]
    assert_array_equal(result, expected)


# This test fails, would be nice to have a proper way to do it.
# By using @xfail, we tell pytest that this test is expected to fail
@pytest.mark.xfail(reason="Reconstruction of edges is not exact.")
def test_get_edges():
    original_a, original_b = load_era5_ab()

    all_levels = list(range(1, 138))
    reconstructed_a, reconstructed_b = get_edges(all_levels)
    assert_almost_equal(original_a, reconstructed_a)
    assert_almost_equal(original_b, reconstructed_b)
