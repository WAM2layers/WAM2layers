from functools import partial

import numpy as np


def pad_boundaries(*args, periodic=False):
    """Add boundary padding to all input arrays.

    Arguments:
        *args: Input arrays to which boundary padding will be added.
        periodic: If True, apply periodic boundary conditions when padding in
            the x-direction. If False, pad with zeros in the x and y directions.

    Returns:
        List of input arrays with boundary padding added.
    """
    zero_xy = partial(np.pad, pad_width=1, mode="constant", constant_values=0)
    periodic_x = partial(np.pad, pad_width=((0, 0), (1, 1)), mode="wrap")
    zero_y = partial(
        np.pad, pad_width=((1, 1), (0, 0)), mode="constant", constant_values=0
    )

    if periodic:
        return [periodic_x(zero_y(arg)) for arg in args]
    return [zero_xy(arg) for arg in args]


def horizontal_advection(s, u, v, periodic_x=False) -> np.ndarray:
    """Calculate advection on a staggered grid using a donor cell scheme.

    Boundaries are padded with 0 or periodic boundaries to maintain the original
    grid size.

    The advection equation reads `grad(u*s)` where u is the velocity vector and
    s is any scalar. Using i for grid cell indices in the x direction and j for
    the y-direction (increasing northwards), the donor cell scheme reads:

    s(i-1/2) = s(i-1) if u(i-1/2) > 0
               s(i+1) otherwise
    s(j-1/2) = s(j-1) if v(j-1/2) > 0
               s(j+1) otherwise

    d(us)/dx = u(i-1/2)*s(i-1/2) - u(i+1/2)*s(i+1/2)
    d(vs)/dy = v(j-1/2)*s(j-1/2) - v(j+1/2)*s(j+1/2)

    adv(s) = d(us)/dx + d(vs)/dy

    Arguments:
        s: array of shape [M, N]. It is assumed that the grid dimensions are
            [latitude, longitude] and latitude is in decreasing order.
        u: array of shape = [M-2, N-1]
        v: array of shape = [M-1, N-2]

    Returns:
        array of shape [M, N]

    Examples:

        Create a simple array:
        >>> s = np.zeros((5, 5))
        >>> s[2, 2] = 1
        >>> u = np.ones((3, 4)) * .5
        >>> v = np.ones((4, 3)) * .5
        >>> s
        array([[0., 0., 0., 0., 0.],
               [0., 0., 0., 0., 0.],
               [0., 0., 1., 0., 0.],
               [0., 0., 0., 0., 0.],
               [0., 0., 0., 0., 0.]])

        Calculate the advection for a single time step (forward)
        >>> horizontal_advection(s, u, v)
        array([[ 0. ,  0. ,  0. ,  0. ,  0. ],
               [ 0. ,  0. ,  0.5,  0. ,  0. ],
               [ 0. ,  0. , -1. ,  0.5,  0. ],
               [ 0. ,  0. ,  0. ,  0. ,  0. ],
               [ 0. ,  0. ,  0. ,  0. ,  0. ]])

        Backtracking can be done by reversing the velocities:
        >>> horizontal_advection(s, -u, -v)
        array([[ 0. ,  0. ,  0. ,  0. ,  0. ],
               [ 0. ,  0. ,  0. ,  0. ,  0. ],
               [ 0. ,  0.5, -1. ,  0. ,  0. ],
               [ 0. ,  0. ,  0.5,  0. ,  0. ],
               [ 0. ,  0. ,  0. ,  0. ,  0. ]])

    """
    # Pad boundaries with ghost cells
    qp, up, vp = pad_boundaries(s, u, v, periodic=periodic_x)
    # shapes: [M+2, N+2], [M, N+1], [M+1, N]

    west = np.s_[:-1]
    east = np.s_[1:]
    south = np.s_[1:]
    north = np.s_[:-1]
    inner = np.s_[1:-1]

    # Donor cell upwind scheme (2 directions seperately)
    uq = np.where(up > 0, up * qp[inner, west], up * qp[inner, east])  # [M, N+1]

    vq = np.where(vp > 0, vp * qp[south, inner], vp * qp[north, inner])  # [M, N+2]

    adv_x = uq[:, west] - uq[:, east]  # [M, N]
    adv_y = vq[south, :] - vq[north, :]  # [M, N]

    return adv_x + adv_y  # [M, N]


def vertical_advection(fv, s_lower, s_upper):
    """Calculate 1d upwind advection of vertical flux.

    Upwind advection with fv positive downwards, so:

    fv * s = fv * s_upper if fv > 0
           = fv * s_lower otherwise
    """
    return np.where(fv >= 0, fv * s_upper, fv * s_lower)


def vertical_dispersion(fv, s_lower, s_upper, kvf):
    """Calculate additional vertical mixing due to convective dispersion.

    dispersion = kvf * |Fv| * dS/dz
    """
    return kvf * np.abs(fv) * (s_upper - s_lower)
