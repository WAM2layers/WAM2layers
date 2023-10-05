from wam2layers.tracking.backtrack import advection
import numpy as np

def test_advection():
    q = np.zeros((5, 5))
    q[2, 2] = 1
    u = np.ones((3, 4)) * .5
    v = np.ones((4, 3)) * .5
    tendency = advection(q, u, v)
    expected = np.array([
        [ 0.,  .5,  0.],
        [ 0., -1.,  .5],
        [ 0.,  0.,  0.]
    ])
    np.testing.assert_almost_equal(tendency, expected)




