"""Generic functions useful for preprocessing various input datasets."""

import numpy as np
from scipy.interpolate import interp1d
import xarray as xr


def resample(variable, divt, count_time, method='interp'):
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


def get_stable_fluxes(fa_e, fa_n, w):
    """Stabilize the outfluxes / influxes.

    During the reduced timestep the water cannot move further than 1/x * the
    gridcell, In other words at least x * the reduced timestep is needed to
    cross a gridcell.
    """
    fa_e_abs = np.abs(fa_e)
    fa_n_abs = np.abs(fa_n)

    stab = 1.0 / 2.0

    # TODO: Check; I don't really understand what's happening here
    fa_e_corrected = (fa_e_abs / (fa_e_abs + fa_n_abs)) * stab * w[:-1, :, :]
    fa_e_stable = np.minimum(fa_e_abs, fa_e_corrected)

    fa_n_corrected = (fa_n_abs / (fa_e_abs + fa_n_abs)) * stab * w[:-1, :, :]
    fa_n_stable = np.minimum(fa_n_abs, fa_n_corrected)

    # get rid of any nan values
    fa_e_stable[np.isnan(fa_e_stable)] = 0
    fa_n_stable[np.isnan(fa_n_stable)] = 0

    # redefine
    fa_e = np.sign(fa_e) * fa_e_stable
    fa_n = np.sign(fa_n) * fa_n_stable

    return fa_e, fa_n


def divergence_zonal(fa_e, periodic_boundary):
    """Define the horizontal fluxes over the zonal boundaries."""
    flux = np.asarray(fa_e)

    fa_e_boundary = np.zeros_like(flux)
    fa_e_boundary[:, :, :-1] = 0.5 * (flux[:, :, :-1] + flux[:, :, 1:])
    if periodic_boundary == True:
        fa_e_boundary[:, :, -1] = 0.5 * (flux[:, :, -1] + flux[:, :, 0])

    fa_w_boundary = np.roll(fa_e_boundary, 1)
    return fa_w_boundary - fa_e_boundary


def divergence_meridional(fa_n):
    """Define the horizontal fluxes over the meridional boundaries."""
    flux = np.asarray(fa_n)

    fa_n_boundary = np.zeros_like(flux)
    fa_n_boundary[:, 1:, :] = 0.5 * (flux[:, :-1, :] + flux[:, 1:, :])

    fa_s_boundary = np.roll(fa_n_boundary, -1, axis=1)
    return fa_s_boundary - fa_n_boundary


def get_vertical_transport(
    fa_e_top,
    fa_e_down,
    fa_n_top,
    fa_n_down,
    e,
    p,
    w_top,
    w_down,
    periodic_boundary
):

    # total moisture in the column
    w = w_top + w_down

    zonal_divergence_top = divergence_zonal(fa_e_top, periodic_boundary)
    zonal_divergence_down = divergence_zonal(fa_e_down, periodic_boundary)
    meridional_divergence_top = divergence_meridional(fa_n_top)
    meridional_divergence_down = divergence_meridional(fa_n_down)

    # check the water balance
    residual_down = np.zeros_like(p)  # residual factor [m3]
    residual_top = np.zeros_like(p)  # residual factor [m3]

    tendency_down = (  # todo why skip the n/s but not w/e edges?
        +zonal_divergence_down[:, 1:-1, :]
        + meridional_divergence_down[:, 1:-1, :]
        - p[:, 1:-1, :] * (w_down[:-1, 1:-1, :] / w[:-1, 1:-1, :])
        + e[:, 1:-1, :]
    )

    tendency_top = (  # todo why skip the n/s but not w/e edges?
        +zonal_divergence_top[:, 1:-1, :]
        + meridional_divergence_top[:, 1:-1, :]
        - p[:, 1:-1, :] * (w_top[:-1, 1:-1, :] / w[:-1, 1:-1, :])
    )

    residual_down[:, 1:-1, :] = (
        w_down[1:, 1:-1, :] - w_down[:-1, 1:-1, :] - tendency_down
    )
    residual_top[:, 1:-1, :] = w_top[1:, 1:-1, :] - w_top[:-1, 1:-1, :] - tendency_top

    # compute the resulting vertical moisture flux; the vertical velocity so
    # that the new residual_down/w_down = residual_top/w_top (positive downward)
    fa_vert_raw = (
        w_down[1:, :, :] / w[1:, :, :] * (residual_down + residual_top) - residual_down
    )

    # stabilize the outfluxes / influxes; during the reduced timestep the
    # vertical flux can maximally empty/fill 1/x of the top or down storage
    stab = 1.0 / 4.0
    fa_vert_stable = np.minimum(
        np.abs(fa_vert_raw), np.minimum(stab * w_top[1:, :, :], stab * w_down[1:, :, :])
    )

    # redefine the vertical flux
    fa_vert = np.sign(fa_vert_raw) * fa_vert_stable

    return fa_vert


def get_grid_info(latitude, longitude):
    """Return grid cell area and lenght sides."""
    dg = 111089.56  # [m] length of 1 degree latitude
    erad = 6.371e6  # [m] Earth radius

    gridcell = np.abs(longitude[1] - longitude[0])  # [degrees] grid cell size

    # new area size calculation:
    lat_n_bound = np.minimum(90.0, latitude + 0.5 * gridcell)
    lat_s_bound = np.maximum(-90.0, latitude - 0.5 * gridcell)

    # TODO check this calculation!
    a_gridcell = np.zeros([len(latitude), 1])
    a_gridcell[:, 0] = (
        (np.pi / 180.0)
        * erad ** 2
        * abs(np.sin(lat_s_bound * np.pi / 180.0) - np.sin(lat_n_bound * np.pi / 180.0))
        * gridcell
    )

    l_ew_gridcell = gridcell * dg  # [m] length eastern/western boundary of a cell
    l_n_gridcell = (
        dg * gridcell * np.cos((latitude + gridcell / 2) * np.pi / 180)
    )  # [m] length northern boundary of a cell
    l_s_gridcell = (
        dg * gridcell * np.cos((latitude - gridcell / 2) * np.pi / 180)
    )  # [m] length southern boundary of a cell
    l_mid_gridcell = 0.5 * (l_n_gridcell + l_s_gridcell)
    return a_gridcell, l_ew_gridcell, l_mid_gridcell


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
