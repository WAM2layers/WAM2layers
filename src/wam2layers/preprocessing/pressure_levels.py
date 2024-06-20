"""Routines specific for preprocessing pressure level data."""
import numpy as np
import xarray as xr

from wam2layers.config import Config
from wam2layers.preprocessing.utils import interpolate, sortby_ndarray


def extend_pressurelevels(
    level_data: xr.Dataset,
    surface_data: xr.Dataset,
    config: Config,
) -> tuple[xr.Dataset, xr.DataArray]:
    """Extend the levels to top of atmosphere, surface and 2-layer boundary level.

    The input data does not have data at the surface and
    and the top of atmosphere. The input data is extended to
    these edges.

    For proper integration, extra layers have to be inserted
    near the 2-layer boundary. This boundary level is hard
    coded as a function of surface pressure and not customizable.

    Args:
        level_data: Dataset with the following variables:
            q: Specific humidity at pressure levels
            u: Eastward horizontal wind speed at pressure levels
            v: Northward horizontal wind speed at pressure levels
        surface_data: Dataset with the following variables:
            ps: Air pressure at the surface
            qs: Specific humidity at the surface
            us: Eastward horizontal wind speed at the surface
            vs: Northward horizontal wind speed at the surface

    Returns:
        Dataset with the following variables (interpolated to new pressure levels):
            p, q, u, v
        The pressure at the division between the two layers
    """
    level_ax: int = level_data["u"].get_axis_num("level")  # type: ignore

    da_p = level_data["level"].broadcast_like(level_data["u"])

    p = da_p.to_numpy()
    u = level_data["u"].to_numpy()
    v = level_data["v"].to_numpy()
    q = level_data["q"].to_numpy()

    # Insert top of atmosphere values
    # Assume wind at top same as values at lowest pressure, humidity at top = 0
    p0 = p.argmin(axis=level_ax).flatten()[0]
    p = np.insert(p, 0, 0.0, axis=level_ax)
    q = np.insert(q, 0, 0.0, axis=level_ax)
    v = np.insert(v, 0, np.take(v, p0, axis=level_ax), axis=level_ax)
    u = np.insert(u, 0, np.take(u, p0, axis=level_ax), axis=level_ax)

    # Insert surface level values
    p = np.insert(p, 0, surface_data["ps"], axis=level_ax)
    q = np.insert(q, 0, surface_data["qs"], axis=level_ax)
    v = np.insert(v, 0, surface_data["vs"], axis=level_ax)
    u = np.insert(u, 0, surface_data["us"], axis=level_ax)

    # Sort arrays by pressure (ascending)
    u = sortby_ndarray(u, p, axis=level_ax)
    v = sortby_ndarray(v, p, axis=level_ax)
    q = sortby_ndarray(q, p, axis=level_ax)
    p = sortby_ndarray(p, p, axis=level_ax)

    # Insert boundary level values (at a dummy pressure value)
    p_boundary = (
        config.pressure_boundary_factor * np.array(surface_data["ps"])
        + config.pressure_boundary_offset
    )
    u = np.insert(u, 0, interpolate(p_boundary, p, u, axis=level_ax), axis=level_ax)
    v = np.insert(v, 0, interpolate(p_boundary, p, v, axis=level_ax), axis=level_ax)
    q = np.insert(q, 0, interpolate(p_boundary, p, q, axis=level_ax), axis=level_ax)
    p = np.insert(p, 0, p_boundary, axis=level_ax)

    # Sort arrays by pressure once more (ascending)
    u = sortby_ndarray(u, p, axis=level_ax)
    v = sortby_ndarray(v, p, axis=level_ax)
    q = sortby_ndarray(q, p, axis=level_ax)
    p = sortby_ndarray(p, p, axis=level_ax)

    # Reset level coordinate as its values have become meaningless
    nlev = np.size(u, axis=level_ax)
    levs = np.arange(nlev)

    # reconstruct dataarrays
    coords = dict(level_data["u"].coords)
    coords["level"] = levs

    output = xr.Dataset()  # New dataset as the coordinate size has changed
    output["u"] = xr.DataArray(u, coords, level_data["u"].dims)
    output["v"] = xr.DataArray(v, coords, level_data["v"].dims)
    output["q"] = xr.DataArray(q, coords, level_data["q"].dims)
    output["p"] = xr.DataArray(p, coords, level_data["u"].dims)
    pb = xr.DataArray(p_boundary, surface_data["ps"].coords, surface_data["ps"].dims)

    return output, pb


def interp_dp_midpoints(
    level_data: xr.Dataset,
) -> tuple[xr.Dataset, xr.DataArray]:
    """Interpolate the data to midpoints to allow for integration to two layers.

    Args:
        level_data: Dataset with the following variables:
            q: Specific humidity at levels
            u: Eastward horizontal wind speed at levels
            v: Northward horizontal wind speed at levels
            p: Air pressure at levels

    Returns:
        Dataset with the following variables:
            p interpolated to midpoints of pressure levels,
            q interpolated to midpoints of pressure levels,
            u interpolated to midpoints of pressure levels,
            v interpolated to midpoints of pressure levels,
        Pressure difference between the pressure level bounds,
    """

    # Calculate pressure jump
    dp = level_data["p"].diff("level")
    assert np.all(dp >= 0), "Pressure levels should increase monotonically"

    # Interpolate to midpoints
    midpoints = 0.5 * (
        level_data["level"].to_numpy()[1:] + level_data["level"].to_numpy()[:-1]
    )
    dp = dp.assign_coords(level=midpoints)
    interped_data = level_data.interp(level=midpoints)
    return interped_data, dp


def mask_below_surface_data(
    level_data: xr.Dataset,
    ps: xr.DataArray,
) -> tuple[xr.Dataset, xr.DataArray]:
    """Mask any points that fall below the surface.

    Args:
        level_data: Dataset with the following variables:
            q: Specific humidity at levels
            u: Eastward horizontal wind speed at levels
            v: Northward horizontal wind speed at levels
            p: Air pressure at levels
        ps: Air pressure at the surface

    Returns:
        The input data, with the values below the surface (based on surface
            air pressure) masked out.
    """

    # mask values below surface
    above_surface = level_data["p"] < ps.broadcast_like(level_data["p"])
    return level_data.where(above_surface)
