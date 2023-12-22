import numpy as np
import pytest
import xarray as xr

from wam2layers.config import BoundingBox
from wam2layers.tracking.io import select_subdomain


@pytest.fixture
def ds():
    lon = np.arange(-180, 180, 30)
    lat = np.arange(90, -90, -30)
    time = np.arange(3)

    p = xr.DataArray(
        np.random.randn(len(time), len(lat), len(lon)),
        coords={"time": time, "latitude": lat, "longitude": lon},
    )
    e = xr.DataArray(
        np.random.randn(len(time), len(lat), len(lon)),
        coords={"time": time, "latitude": lat, "longitude": lon},
    )
    return xr.Dataset({"e": e, "p": p})


def test_select_subdomain(ds):
    """Check subdomain that fits in [-180, 180]."""
    bbox = BoundingBox(0, -20, 90, 60)
    sub = select_subdomain(ds, bbox)
    assert np.array_equal(sub.latitude, [60, 30, 0])
    assert np.array_equal(sub.longitude, [0, 30, 60, 90])


def test_rolling_subdomain(ds):
    """Check subdomain that does not fit in [-180, 180] or [0, 360]."""
    bbox = BoundingBox(130, -20, 40, 60)
    sub = select_subdomain(ds, bbox)
    assert np.array_equal(sub.longitude, [150, -180, -150, -120, -90, -60, -30, 0, 30])
