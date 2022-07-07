"""Generic functions useful for preprocessing various input datasets."""

import numpy as np
from scipy.interpolate import interp1d
import xarray as xr


def resample_ec_earth(variable, divt, count_time, method='interp'):
    """Resample the variable to a given number of timesteps."""

    _, nlat, nlon = variable.shape
    ntime = count_time * divt
    shape = (ntime, nlat, nlon)
    new_var = np.nan * np.zeros(shape)

    oddvector = np.zeros(ntime)
    partvector = np.zeros(ntime)
    for o in np.arange(0, ntime, divt):
        for i in range(1, divt):
            oddvector[o + i - 1] = (divt - i) / divt
            partvector[o + i] = i / divt

    for t in range(ntime):
        idx = int((t+1) / divt + oddvector[t])

        if method == 'bfill':
            new_var[t] = (1 / divt) * variable[idx - 1]
        elif method == 'interp':
            new_var[t] = variable[idx - 1] + partvector[t] * (variable[idx] - variable[idx - 1])
        else:
            raise ValueError(f"Unknown resample method {method}")

    return new_var


def stabilize_fluxes(fx, fy, s):
    """Stabilize the outfluxes / influxes.

    During the reduced timestep the water cannot move further than 1/x * the
    gridcell, In other words at least x * the reduced timestep is needed to
    cross a gridcell.
    """
    fx_abs = np.abs(fx)
    fy_abs = np.abs(fy)
    ft_abs = fx_abs + fy_abs

    # TODO: Check; I don't really understand what's happening here
    fx_corrected = 1/2 * fx_abs / ft_abs * s[:-1, :, :]
    fx_stable = np.minimum(fx_abs, fx_corrected)

    fy_corrected = 1/2 * fy_abs / ft_abs * s[:-1, :, :]
    fy_stable = np.minimum(fy_abs, fy_corrected)

    # get rid of any nan values
    fx_stable[np.isnan(fx_stable)] = 0
    fy_stable[np.isnan(fy_stable)] = 0

    # redefine
    fx = np.sign(fx) * fx_stable
    fy = np.sign(fy) * fy_stable

    return fx, fy


def divergence_zonal(fx, periodic_boundary):
    """Define the horizontal fluxes over the zonal boundaries."""
    flux = np.asarray(fx)

    fe_boundary = np.zeros_like(flux)
    fe_boundary[:, :, :-1] = 0.5 * (flux[:, :, :-1] + flux[:, :, 1:])
    if periodic_boundary:
        fe_boundary[:, :, -1] = 0.5 * (flux[:, :, -1] + flux[:, :, 0])

    fw_boundary = np.roll(fe_boundary, 1)
    return fw_boundary - fe_boundary


def divergence_meridional(fy):
    """Define the horizontal fluxes over the meridional boundaries."""
    flux = np.asarray(fy)

    fn_boundary = np.zeros_like(flux)
    fn_boundary[:, 1:, :] = 0.5 * (flux[:, :-1, :] + flux[:, 1:, :])

    fs_boundary = np.roll(fn_boundary, -1, axis=1)
    return fs_boundary - fn_boundary


def get_vertical_transport(
    fx_upper,
    fx_lower,
    fy_upper,
    fy_lower,
    evap,
    precip,
    w_upper,
    w_lower,
    periodic_boundary,
    kvf
):

    # total moisture in the column
    w = w_upper + w_lower

    zonal_divergence_upper = divergence_zonal(fx_upper, periodic_boundary)
    zonal_divergence_lower = divergence_zonal(fx_lower, periodic_boundary)
    meridional_divergence_upper = divergence_meridional(fy_upper)
    meridional_divergence_lower = divergence_meridional(fy_lower)

    # check the water balance
    residual_lower = np.zeros_like(precip)  # residual factor [m3]
    residual_upper = np.zeros_like(precip)  # residual factor [m3]

    tendency_lower = (  # todo why skip the n/s but not w/e edges?
        +zonal_divergence_lower[:, 1:-1, :]
        + meridional_divergence_lower[:, 1:-1, :]
        - precip[:, 1:-1, :] * (w_lower[:-1, 1:-1, :] / w[:-1, 1:-1, :])
        + evap[:, 1:-1, :]
    )

    tendency_upper = (  # todo why skip the n/s but not w/e edges?
        +zonal_divergence_upper[:, 1:-1, :]
        + meridional_divergence_upper[:, 1:-1, :]
        - precip[:, 1:-1, :] * (w_upper[:-1, 1:-1, :] / w[:-1, 1:-1, :])
    )

    residual_lower[:, 1:-1, :] = (
        w_lower[1:, 1:-1, :] - w_lower[:-1, 1:-1, :] - tendency_lower
    )
    residual_upper[:, 1:-1, :] = w_upper[1:, 1:-1, :] - w_upper[:-1, 1:-1, :] - tendency_upper

    # compute the resulting vertical moisture flux; the vertical velocity so
    # that the new residual_lower/w_lower = residual_upper/w_upper (positive downward)
    fa_vert_raw = (
        w_lower[1:, :, :] / w[1:, :, :] * (residual_lower + residual_upper) - residual_lower
    )

    # stabilize the outfluxes / influxes; during the reduced timestep the
    # vertical flux can maximally empty/fill 1/x of the top or down storage
    stab = 1.0 / (kvf + 1.0)
    fa_vert_stable = np.minimum(
        np.abs(fa_vert_raw), np.minimum(stab * w_upper[1:, :, :], stab * w_lower[1:, :, :])
    )

    # redefine the vertical flux
    fa_vert = np.sign(fa_vert_raw) * fa_vert_stable

    return fa_vert


def get_grid_info(ds):
    """Return grid cell area and lenght of the sides."""
    dg = 111089.56  # [m] length of 1 degree latitude
    erad = 6.371e6  # [m] Earth radius

    latitude = ds.latitude.values
    longitude = ds.longitude.values
    grid_spacing = np.abs(longitude[1] - longitude[0])  # [degrees]

    # Calculate area TODO check this calculation!
    lat_n = np.minimum(90.0, latitude + 0.5 * grid_spacing)
    lat_s = np.maximum(-90.0, latitude - 0.5 * grid_spacing)

    a = np.pi / 180.0 * erad ** 2 * grid_spacing * abs(
        np.sin(lat_s * np.pi / 180.0) - np.sin(lat_n * np.pi / 180.0)
    )

    # Calculate faces
    ly = grid_spacing * dg  # [m] length eastern/western boundary of a cell
    lx_n_gridcell = ly * np.cos((latitude + grid_spacing / 2) * np.pi / 180)
    lx_s_gridcell = ly * np.cos((latitude - grid_spacing / 2) * np.pi / 180)
    lx = 0.5 * (lx_n_gridcell + lx_s_gridcell)
    return a, ly, lx


def join_levels(pressure_level_data, surface_level_data):
    """Combine 3d pressure level and 2d surface level data.

    A dummy value of 1000 hPa is inserted for the surface pressure
    """
    bottom = pressure_level_data.isel(lev=0).copy()
    bottom.values = surface_level_data.values
    bottom["lev"] = 100000.0  # dummy value

    # NB: plevs needs to come first to preserve dimension order
    return xr.concat([pressure_level_data, bottom], dim="lev").sortby(
        "lev", ascending=False
    )


def repeat_upper_level(pressure_level_data, fill_value=None):
    """Add one more level with the same value as the upper level.

    A dummy value of 100 hPa is inserted for the top level pressure
    """
    top_level_data = pressure_level_data.isel(lev=-1).copy()
    top_level_data["lev"] = 10000.0

    if fill_value is not None:
        top_level_data.values = np.ones_like(top_level_data) * fill_value

    return xr.concat([pressure_level_data, top_level_data], dim="lev")


def get_new_target_levels(surface_pressure, p_boundary, n_levels=40):
    """Build a numpy array with new pressure levels including the boundary."""
    # remove xarray labels if present
    surface_pressure = np.array(surface_pressure)

    pressure_upper = 20000
    dp = (surface_pressure - pressure_upper) / (n_levels - 1)

    levels = np.arange(0, n_levels)[None, :, None, None]
    ntime, _, nlat, nlon = surface_pressure.shape

    # Note the extra layer of zeros at the bottom and top of the array
    # TODO: find out why
    new_p = np.zeros((ntime, n_levels + 2, nlat, nlon))
    new_p[:, 1:-1, :, :] = surface_pressure - dp * levels

    mask = np.where(new_p > p_boundary, 1, 0)
    mask[:, 0, :, :] = 1  # bottom value is always 1

    # Insert the boundary in the new levels, pushing higher layers one index up
    # e.g. with the boundary at 850 hPa:
    # [0, 1000, 900, 800, 700, 600, 0] --> [0, 1000, 900, 850, 800, 700, 600]
    new_p[:, :-1, :, :] = (
        mask[:, 1:, :, :] * new_p[:, 1:, :, :]
        + (1 - mask[:, 1:, :, :]) * new_p[:, :-1, :, :]
    )

    new_p[:, 1:, :, :] = np.where(
        new_p[:, :-1, :, :] == new_p[:, 1:, :, :],
        p_boundary,
        new_p[:, 1:, :, :],
    )
    return new_p


def interpolate(old_var, old_pressure_levels, new_pressure_levels, type="linear"):
    """Interpolate old_var to new_pressure_levels."""
    new_var = np.zeros_like(new_pressure_levels)

    ntime, _, nlat, nlon = old_var.shape

    for t in range(ntime):
        for i in range(nlat):
            for j in range(nlon):

                pressure_1d = old_pressure_levels[t, :, i, j]
                var_1d = old_var[t, :, i, j]

                pressure_1d = pressure_1d[~pressure_1d.mask]
                var_1d = var_1d[~var_1d.mask]

                f_q = interp1d(pressure_1d, var_1d, type)
                new_var[t, :, i, j] = f_q(new_pressure_levels[t, :, i, j])

    return new_var
