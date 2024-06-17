"""Functionality related to the grid."""

import numpy as np
import xarray as xr


def stagger_x(f):
    """Interpolate f to the grid cell interfaces.
    Only the values at the interior interfaces are returned
    Arguments:
        f: 2d array of shape [M, N]
    Returns:
        2d array of shape [M-2, N-1]
    """
    return 0.5 * (f[:, :-1] + f[:, 1:])[1:-1, :]


def stagger_y(f):
    """Interpolate f to the grid cell interfaces.
    Only the values at the interior interfaces are returned
    Arguments:
        f: 2d array of shape [M, N]
    Returns:
        2d array of shape [M-1, N-2]
    """
    return 0.5 * (f[:-1, :] + f[1:, :])[:, 1:-1]


def get_grid_info(ds):
    """Return grid cell area and lenght of the sides."""
    dg = 111089.56  # [m] length of 1 degree latitude
    erad = 6.371e6  # [m] Earth radius

    latitude = ds.latitude.values
    longitude = ds.longitude.values
    dx_deg = np.abs(longitude[1] - longitude[0])  # [degrees]
    dy_deg = np.abs(latitude[1] - latitude[0])  # [degrees]

    # Calculate area TODO check this calculation!
    lat_n = np.minimum(90.0, latitude + 0.5 * dx_deg)
    lat_s = np.maximum(-90.0, latitude - 0.5 * dx_deg)

    a = (
        np.pi
        / 180.0
        * erad**2
        * dx_deg
        * abs(np.sin(lat_s * np.pi / 180.0) - np.sin(lat_n * np.pi / 180.0))
    )

    # Calculate faces
    dy = dy_deg * dg  # [m] grid spacing in meridional direction
    dx_n_gridcell = dx_deg * dg * np.cos((latitude + dx_deg / 2) * np.pi / 180)
    dx_s_gridcell = dx_deg * dg * np.cos((latitude - dx_deg / 2) * np.pi / 180)
    dx = 0.5 * (dx_n_gridcell + dx_s_gridcell)  # [m] grid spacing in zonal direction
    return a, dy, dx[:, None]


def get_boundary(field, periodic=False):
    """Return a mask with 1 along the boundary and 0 in the interior."""
    boundary = xr.ones_like(field, dtype=bool)
    # Mask interior
    if periodic:
        boundary[1:-1, :] = 0
    else:
        boundary[1:-1, 1:-1] = 0

    return boundary
