from functools import partial

import numpy as np

from wam2layers.config import Config
from wam2layers.utils.profiling import ProgressTracker


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

    vq = np.where(vp > 0, vp * qp[south, inner], vp * qp[north, inner])  # [M+1, N]

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


def stabilize_fluxes(F, S, progress_tracker: ProgressTracker, config: Config, t):
    """Stabilize the outfluxes / influxes.

    CFL: Water cannot move further than one grid cell per timestep.

    Arguments:
        F: xr.Dataset holding the fluxes
        S: xr.Dataset holding the states
        progress_tracker: ProgressTracker object for writing log messages
        config: Config object with the configuration for this experiment
        t: current time (for writing output)
    """
    for level in ["upper", "lower"]:
        fx = F["fx_" + level] * config.timestep
        fy = F["fy_" + level] * config.timestep
        s = S["s_" + level]

        fx_abs = np.abs(fx)
        fy_abs = np.abs(fy)
        ft_abs = fx_abs + fy_abs

        # TODO: make 1/2 configurable?
        fx_limit = 1 / 2 * fx_abs / ft_abs * s.values
        fx_stable = np.minimum(fx_abs, fx_limit)

        fy_limit = 1 / 2 * fy_abs / ft_abs * s.values
        fy_stable = np.minimum(fy_abs, fy_limit)

        progress_tracker.track_stability_correction(fy_stable, fy_abs, config, t)

        # Get rid of any nan values
        fx_stable.fillna(0)
        fy_stable.fillna(0)

        # Re-instate the sign and convert back to flux instead of accumulation
        F["fx_" + level] = np.sign(fx) * fx_stable / config.timestep
        F["fy_" + level] = np.sign(fy) * fy_stable / config.timestep


def divergence(fx, fy):
    # Note: latitude decreasing, hence negative fy gradient is divergence
    return np.gradient(fx, axis=-1) - np.gradient(fy, axis=-2)


def calculate_fz(F, S0, S1, dt, kvf):
    """Calculate the vertical fluxes.

    The vertical flux is calculated as a closure term. Residuals are distributed
    proportionally to the amount of moisture in each layer.

    The flux is constrained such that it can never exceed

    Arguments:
        F: xarray dataset with fluxes evaluated at temporal midpoints between states
        S0: xarray dataset with states at current time t
        S1: xarray dataset with states at updated time t+1 (always forward looking)
        dt: timestep
        kvf: net to gross vertical flux multiplication parameter. With the a value of 0,
            the net and gross fluxes are equal.

    Returns:
        fz: vertical flux, positive downward

    Examples:

        Create dummy input data. In the absence of any horizontal fluxes and
        sources or sinks, this result can be verified manually.

        >>> import numpy as np
        >>> import xarray as xr
        >>> from wam2layers.preprocessing.shared import calculate_fz
        >>>
        >>> F = xr.Dataset({
        ...     'precip': xr.DataArray(np.zeros(((3, 3))), dims=['lat', 'lon']),
        ...     'evap': xr.DataArray(np.zeros(((3, 3))), dims=['lat', 'lon']),
        ...     'fx_upper': xr.DataArray(np.zeros(((3, 3))), dims=['lat', 'lon']),
        ...     'fy_upper': xr.DataArray(np.zeros(((3, 3))), dims=['lat', 'lon']),
        ...     'fx_lower': xr.DataArray(np.zeros(((3, 3))), dims=['lat', 'lon']),
        ...     'fy_lower': xr.DataArray(np.zeros(((3, 3))), dims=['lat', 'lon']),
        ... })
        >>> S0 = xr.Dataset({
        ...     's_upper': 10 * xr.DataArray(np.ones(((3, 3))), dims=['lat', 'lon']),
        ...     's_lower': 6 * xr.DataArray(np.ones(((3, 3))), dims=['lat', 'lon']),
        ... })
        >>> S1 = xr.Dataset({
        ...     's_upper': 8 * xr.DataArray(np.ones(((3, 3))), dims=['lat', 'lon']),
        ...     's_lower': 7 * xr.DataArray(np.ones(((3, 3))), dims=['lat', 'lon']),
        ... })
        >>> calculate_fz(F, S0, S1, dt=1)
        <xarray.DataArray (lat: 3, lon: 3)>
        array([[1.46666667, 1.46666667, 1.46666667],
               [1.46666667, 1.46666667, 1.46666667],
               [1.46666667, 1.46666667, 1.46666667]])
        Dimensions without coordinates: lat, lon

    """
    s_mean = (S1 + S0) / 2
    s_total = s_mean.s_upper + s_mean.s_lower
    s_rel = s_mean / s_total

    # Evaluate all terms in the moisture balance execpt the unknown Fz and err
    foo_upper = (S1 - S0).s_upper + dt * (
        divergence(F.fx_upper, F.fy_upper) + F.precip.values * s_rel.s_upper
    )
    foo_lower = (S1 - S0).s_lower + dt * (
        divergence(F.fx_lower, F.fy_lower) + F.precip.values * s_rel.s_lower - F.evap
    )

    # compute the resulting vertical moisture flux; the vertical velocity so
    # that the new err_lower/s_lower = err_upper/s_upper (positive downward)
    fz = foo_lower - S1.s_lower / (S1.s_lower + S1.s_upper) * (foo_upper + foo_lower)

    # TODO: verify that err_lower/s_lower = err_upper/s_upper is satisfied on debug log

    # stabilize the outfluxes / influxes; during the reduced timestep the
    # vertical flux can maximally empty/fill 1/x of the top or down storage
    stab = 1.0 / (kvf + 1.0)
    flux_limit = np.minimum(
        s_mean.s_upper, s_mean.s_lower
    )  # TODO why is this not 'just' the upstream bucket?
    fz_stable = np.minimum(np.abs(fz), stab * flux_limit)

    # Reinstate the sign and convert accumulation back to flux
    return np.sign(fz) * fz_stable / dt
