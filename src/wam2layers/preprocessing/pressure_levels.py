from wam2layers.preprocessing.utils import interpolate, sortby_ndarray


import numpy as np
import xarray as xr


def extend_pressurelevels(
    input: xr.Dataset,
    surface: xr.Dataset,
) -> tuple[xr.Dataset, xr.DataArray]:
    """Extend the levels to top of atmosphere, surface and 2-layer boundary level.

    The input data does not have data at the surface and
    and the top of atmosphere. The input data is extended to
    these edges.

    For proper integration, extra layers have to be inserted
    near the 2-layer boundary. This boundary level is hard
    coded as a function of surface pressure and not customizable.

    Args:
        input_data: Dataset with the following variables:
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
    level_ax: int = input["u"].get_axis_num("level")  # type: ignore

    da_p = input["level"].broadcast_like(input["u"]) * 100  # Pa
    da_p.attrs["units"] = "Pa"

    p = da_p.to_numpy()
    u = input["u"].to_numpy()
    v = input["v"].to_numpy()
    q = input["q"].to_numpy()

    # Insert top of atmosphere values
    # Assume wind at top same as values at lowest pressure, humidity at top = 0
    p0 = p.argmin(axis=level_ax).flatten()[
        0
    ]  # index along level axis of minimum pressure
    p = np.insert(p, 0, 0.0, axis=level_ax)
    q = np.insert(q, 0, 0.0, axis=level_ax)
    v = np.insert(v, 0, np.take(v, p0, axis=level_ax), axis=level_ax)
    u = np.insert(u, 0, np.take(u, p0, axis=level_ax), axis=level_ax)

    # Insert surface level values
    p = np.insert(p, 0, surface["ps"], axis=level_ax)
    q = np.insert(q, 0, surface["qs"], axis=level_ax)
    v = np.insert(v, 0, surface["vs"], axis=level_ax)
    u = np.insert(u, 0, surface["us"], axis=level_ax)

    # Sort arrays by pressure (ascending)
    u = sortby_ndarray(u, p, axis=level_ax)
    v = sortby_ndarray(v, p, axis=level_ax)
    q = sortby_ndarray(q, p, axis=level_ax)
    p = sortby_ndarray(p, p, axis=level_ax)

    # Insert boundary level values (at a dummy pressure value)
    p_boundary = 0.72878581 * np.array(surface["ps"]) + 7438.803223
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
    coords = dict(input["u"].coords)
    coords["level"] = levs

    output = xr.Dataset()
    output["u"] = xr.DataArray(u, coords, input["u"].dims, attrs=input["u"].attrs)
    output["v"] = xr.DataArray(v, coords, input["v"].dims, attrs=input["v"].attrs)
    output["q"] = xr.DataArray(q, coords, input["q"].dims, attrs=input["q"].attrs)
    output["p"] = xr.DataArray(p, coords, input["u"].dims, attrs=surface["ps"].attrs)
    pb = xr.DataArray(p_boundary, surface["ps"].coords, surface["ps"].dims)

    return output, pb


def interp_dp_midpoints(
    extended_data: xr.Dataset,
    ps: xr.DataArray,
):
    """Interpolate the data to midpoints to allow for integration to two layers.

    Args:
        Dataset with the following variables:
            q: Specific humidity at levels
            u: Eastward horizontal wind speed at levels
            v: Northward horizontal wind speed at levels
            p: Air pressure at levels
            ps: Air pressure at the surface

    Returns:
        Pressure difference between the pressure level bounds,
        p interpolated to midpoints of pressure levels,
        q interpolated to midpoints of pressure levels,
        u interpolated to midpoints of pressure levels,
        v interpolated to midpoints of pressure levels,
    """

    # Calculate pressure jump
    dp = extended_data["p"].diff("level")
    assert np.all(dp >= 0), "Pressure levels should increase monotonically"

    # Interpolate to midpoints
    midpoints = 0.5 * (
        extended_data["level"].to_numpy()[1:] + 
        extended_data["level"].to_numpy()[:-1]
    )
    dp = dp.assign_coords(level=midpoints)
    interped_data = extended_data.interp(level=midpoints)
    
    # mask values below surface
    above_surface = interped_data["p"] < ps.broadcast_like(interped_data["p"])
    interped_data = interped_data.where(above_surface)

    return interped_data, dp