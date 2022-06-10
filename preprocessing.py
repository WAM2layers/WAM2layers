"""Generic functions useful for preprocessing various input datasets."""

import numpy as np


def getrefined_new(
    fa_e_top,
    fa_n_top,
    fa_e_down,
    fa_n_down,
    w_top,
    w_down,
    e,
    p,
    divt,
    count_time,
    latitude,
    longitude,
):
    """This definition refines the timestep of the data."""
    # within this new definition of refined i do a linear interpolation over time of my fluxes
    # Imme: change the timesteps from 6-hourly and 3-hourly to 96 timesteps a day

    # for 3 hourly information
    divt2 = divt / 2.0
    oddvector2 = np.zeros((1, int(count_time * 2 * divt2)))
    partvector2 = np.zeros((1, int(count_time * 2 * divt2)))
    da = np.arange(1, divt2)

    for o in np.arange(0, int(count_time * 2 * divt2), int(divt2)):
        for i in range(len(da)):
            oddvector2[0, o + i] = (divt2 - da[i]) / divt2
            partvector2[0, o + i + 1] = da[i] / divt2

    e_small = np.nan * np.zeros(
        (int(count_time * 2 * divt2), len(latitude), len(longitude))
    )
    for t in range(1, int(count_time * 2 * divt2) + 1):
        e_small[t - 1] = (1.0 / divt2) * e[int(t / divt2 + oddvector2[0, t - 1] - 1)]
    e = e_small

    p_small = np.nan * np.zeros(
        (int(count_time * 2 * divt2), len(latitude), len(longitude))
    )
    for t in range(1, int(count_time * 2 * divt2) + 1):
        p_small[t - 1] = (1.0 / divt2) * p[int(t / divt2 + oddvector2[0, t - 1] - 1)]
    p = p_small

    # for 6 hourly info
    oddvector = np.zeros((1, int(count_time * divt)))
    partvector = np.zeros((1, int(count_time * divt)))
    da = np.arange(1, divt)
    divt = float(divt)
    for o in np.arange(0, int(count_time * divt), int(divt)):
        for i in range(len(da)):
            oddvector[0, o + i] = (divt - da[i]) / divt
            partvector[0, o + i + 1] = da[i] / divt

    w_top_small = np.nan * np.zeros(
        (int(count_time * divt + 1), len(latitude), len(longitude))
    )
    w_down_small = np.nan * np.zeros(
        (int(count_time * divt + 1), len(latitude), len(longitude))
    )

    fa_e_down_small = np.nan * np.zeros(
        (int(count_time * divt), len(latitude), len(longitude))
    )
    fa_n_down_small = np.nan * np.zeros(
        (int(count_time * divt), len(latitude), len(longitude))
    )
    fa_e_top_small = np.nan * np.zeros(
        (int(count_time * divt), len(latitude), len(longitude))
    )
    fa_n_top_small = np.nan * np.zeros(
        (int(count_time * divt), len(latitude), len(longitude))
    )

    for t in range(1, int(count_time * divt) + 1):
        w_top_small[t - 1] = w_top[
            int(t / divt + oddvector[0, t - 1] - 1)
        ] + partvector[0, t - 1] * (
            w_top[int(t / divt + oddvector[0, t - 1])]
            - w_top[int(t / divt + oddvector[0, t - 1] - 1)]
        )
        w_top_small[-1] = w_top[-1]
        w_down_small[t - 1] = w_down[
            int(t / divt + oddvector[0, t - 1] - 1)
        ] + partvector[0, t - 1] * (
            w_down[int(t / divt + oddvector[0, t - 1])]
            - w_down[int(t / divt + oddvector[0, t - 1] - 1)]
        )
        w_down_small[-1] = w_down[-1]

        fa_e_down_small[t - 1] = fa_e_down[
            int(t / divt + oddvector[0, t - 1] - 1)
        ] + partvector[0, t - 1] * (
            fa_e_down[int(t / divt + oddvector[0, t - 1])]
            - fa_e_down[int(t / divt + oddvector[0, t - 1] - 1)]
        )
        fa_n_down_small[t - 1] = fa_n_down[
            int(t / divt + oddvector[0, t - 1] - 1)
        ] + partvector[0, t - 1] * (
            fa_n_down[int(t / divt + oddvector[0, t - 1])]
            - fa_n_down[int(t / divt + oddvector[0, t - 1] - 1)]
        )
        fa_e_top_small[t - 1] = fa_e_top[
            int(t / divt + oddvector[0, t - 1] - 1)
        ] + partvector[0, t - 1] * (
            fa_e_top[int(t / divt + oddvector[0, t - 1])]
            - fa_e_top[int(t / divt + oddvector[0, t - 1] - 1)]
        )
        fa_n_top_small[t - 1] = fa_n_top[
            int(t / divt + oddvector[0, t - 1] - 1)
        ] + partvector[0, t - 1] * (
            fa_n_top[int(t / divt + oddvector[0, t - 1])]
            - fa_n_top[int(t / divt + oddvector[0, t - 1] - 1)]
        )

    w_top = w_top_small
    w_down = w_down_small
    fa_e_down = fa_e_down_small
    fa_n_down = fa_n_down_small
    fa_e_top = fa_e_top_small
    fa_n_top = fa_n_top_small

    return fa_e_top, fa_n_top, fa_e_down, fa_n_down, e, p, w_top, w_down


def get_stable_fluxes(fa_e, fa_n, w):
    """Stabilize the outfluxes / influxes.

    During the reduced timestep the water cannot move further than 1/x * the
    gridcell, In other words at least x * the reduced timestep is needed to
    cross a gridcell.
    """
    fa_e_abs = np.abs(fa_e)
    fa_n_abs = np.abs(fa_n)

    stab = 1.0 / 2.0

    fa_e_corrected = (fa_e_abs / (fa_e_abs + fa_n_abs)) * stab * w[:-1, :, :]
    fa_e_stable = np.minimum(fa_e_abs, fa_e_corrected)

    fa_n_corrected = (fa_n_abs / (fa_e_abs + fa_n_abs)) * stab * w[:-1, :, :]
    fa_n_stable = np.minimum(fa_n_abs, fa_n_corrected)

    # get rid of the nan values
    fa_e_stable[np.isnan(fa_e_stable)] = 0
    fa_n_stable[np.isnan(fa_n_stable)] = 0

    # redefine
    fa_e = np.sign(fa_e) * fa_e_stable
    fa_n = np.sign(fa_n) * fa_n_stable

    return fa_e, fa_n


def divergence_zonal(fa_e):
    """Define the horizontal fluxes over the zonal boundaries."""
    flux = np.asarray(fa_e)

    fa_e_boundary = np.zeros_like(flux)
    fa_e_boundary[:, :, :-1] = 0.5 * (flux[:, :, :-1] + flux[:, :, 1:])

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
):

    # total moisture in the column
    w = w_top + w_down

    zonal_divergence_top = divergence_zonal(fa_e_top)
    zonal_divergence_down = divergence_zonal(fa_e_down)
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
