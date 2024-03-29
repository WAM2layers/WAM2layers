import numpy as np

from wam2layers.tracking.core import horizontal_advection


def test_advection():
    q = np.zeros((5, 5))
    q[2, 2] = 1
    u = np.ones((3, 4)) * 0.5
    v = np.ones((4, 3)) * 0.5
    tendency = horizontal_advection(q, u, v)
    expected = np.array([[0.0, 0.5, 0.0], [0.0, -1.0, 0.5], [0.0, 0.0, 0.0]])
    np.testing.assert_almost_equal(tendency, expected)
