"""Generic functions useful for preprocessing various input datasets."""

import numpy as np


# within this new definition of refined I do a linear interpolation over time of my fluxes
def getrefined_new(
    Fa_E_top,
    Fa_N_top,
    Fa_E_down,
    Fa_N_down,
    W_top,
    W_down,
    E,
    P,
    divt,
    count_time,
    latitude,
    longitude,
):
    # This definition refines the timestep of the data
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

    E_small = np.nan * np.zeros(
        (int(count_time * 2 * divt2), len(latitude), len(longitude))
    )
    for t in range(1, int(count_time * 2 * divt2) + 1):
        E_small[t - 1] = (1.0 / divt2) * E[int(t / divt2 + oddvector2[0, t - 1] - 1)]
    E = E_small

    P_small = np.nan * np.zeros(
        (int(count_time * 2 * divt2), len(latitude), len(longitude))
    )
    for t in range(1, int(count_time * 2 * divt2) + 1):
        P_small[t - 1] = (1.0 / divt2) * P[int(t / divt2 + oddvector2[0, t - 1] - 1)]
    P = P_small

    # for 6 hourly info
    oddvector = np.zeros((1, int(count_time * divt)))
    partvector = np.zeros((1, int(count_time * divt)))
    da = np.arange(1, divt)
    divt = float(divt)
    for o in np.arange(0, int(count_time * divt), int(divt)):
        for i in range(len(da)):
            oddvector[0, o + i] = (divt - da[i]) / divt
            partvector[0, o + i + 1] = da[i] / divt

    W_top_small = np.nan * np.zeros(
        (int(count_time * divt + 1), len(latitude), len(longitude))
    )
    W_down_small = np.nan * np.zeros(
        (int(count_time * divt + 1), len(latitude), len(longitude))
    )

    Fa_E_down_small = np.nan * np.zeros(
        (int(count_time * divt), len(latitude), len(longitude))
    )
    Fa_N_down_small = np.nan * np.zeros(
        (int(count_time * divt), len(latitude), len(longitude))
    )
    Fa_E_top_small = np.nan * np.zeros(
        (int(count_time * divt), len(latitude), len(longitude))
    )
    Fa_N_top_small = np.nan * np.zeros(
        (int(count_time * divt), len(latitude), len(longitude))
    )

    for t in range(1, int(count_time * divt) + 1):
        W_top_small[t - 1] = W_top[
            int(t / divt + oddvector[0, t - 1] - 1)
        ] + partvector[0, t - 1] * (
            W_top[int(t / divt + oddvector[0, t - 1])]
            - W_top[int(t / divt + oddvector[0, t - 1] - 1)]
        )
        W_top_small[-1] = W_top[-1]
        W_down_small[t - 1] = W_down[
            int(t / divt + oddvector[0, t - 1] - 1)
        ] + partvector[0, t - 1] * (
            W_down[int(t / divt + oddvector[0, t - 1])]
            - W_down[int(t / divt + oddvector[0, t - 1] - 1)]
        )
        W_down_small[-1] = W_down[-1]

        Fa_E_down_small[t - 1] = Fa_E_down[
            int(t / divt + oddvector[0, t - 1] - 1)
        ] + partvector[0, t - 1] * (
            Fa_E_down[int(t / divt + oddvector[0, t - 1])]
            - Fa_E_down[int(t / divt + oddvector[0, t - 1] - 1)]
        )
        Fa_N_down_small[t - 1] = Fa_N_down[
            int(t / divt + oddvector[0, t - 1] - 1)
        ] + partvector[0, t - 1] * (
            Fa_N_down[int(t / divt + oddvector[0, t - 1])]
            - Fa_N_down[int(t / divt + oddvector[0, t - 1] - 1)]
        )
        Fa_E_top_small[t - 1] = Fa_E_top[
            int(t / divt + oddvector[0, t - 1] - 1)
        ] + partvector[0, t - 1] * (
            Fa_E_top[int(t / divt + oddvector[0, t - 1])]
            - Fa_E_top[int(t / divt + oddvector[0, t - 1] - 1)]
        )
        Fa_N_top_small[t - 1] = Fa_N_top[
            int(t / divt + oddvector[0, t - 1] - 1)
        ] + partvector[0, t - 1] * (
            Fa_N_top[int(t / divt + oddvector[0, t - 1])]
            - Fa_N_top[int(t / divt + oddvector[0, t - 1] - 1)]
        )

    W_top = W_top_small
    W_down = W_down_small
    Fa_E_down = Fa_E_down_small
    Fa_N_down = Fa_N_down_small
    Fa_E_top = Fa_E_top_small
    Fa_N_top = Fa_N_top_small

    return Fa_E_top, Fa_N_top, Fa_E_down, Fa_N_down, E, P, W_top, W_down


# Code
def change_units(
    Fa_E_top_kgpmps,
    Fa_E_down_kgpmps,
    Fa_N_top_kgpmps,
    Fa_N_down_kgpmps,
    timestep,
    divt,
    L_EW_gridcell,
    density_water,
    L_N_gridcell,
    L_S_gridcell,
    latitude,
):

    # convert to m3
    Fa_E_top_m3 = (
        Fa_E_top_kgpmps * timestep / float(divt) * L_EW_gridcell / density_water
    )  # [kg*m^-1*s^-1*s*m*kg^-1*m^3]=[m3]
    Fa_E_down_m3 = (
        Fa_E_down_kgpmps * timestep / float(divt) * L_EW_gridcell / density_water
    )  # [s*m*kg*m^-1*s^-1*kg^-1*m^3]=[m3]

    Fa_N_top_swap = np.zeros(
        (len(latitude), int(count_time * float(divt)), len(longitude))
    )
    Fa_N_down_swap = np.zeros(
        (len(latitude), int(count_time * float(divt)), len(longitude))
    )
    Fa_N_top_kgpmps_swap = np.swapaxes(Fa_N_top_kgpmps, 0, 1)
    Fa_N_down_kgpmps_swap = np.swapaxes(Fa_N_down_kgpmps, 0, 1)
    for c in range(len(latitude)):
        Fa_N_top_swap[c] = (
            Fa_N_top_kgpmps_swap[c]
            * timestep
            / float(divt)
            * 0.5
            * (L_N_gridcell[c] + L_S_gridcell[c])
            / density_water
        )  # [s*m*kg*m^-1*s^-1*kg^-1*m^3]=[m3]
        Fa_N_down_swap[c] = (
            Fa_N_down_kgpmps_swap[c]
            * timestep
            / float(divt)
            * 0.5
            * (L_N_gridcell[c] + L_S_gridcell[c])
            / density_water
        )  # [s*m*kg*m^-1*s^-1*kg^-1*m^3]=[m3]

    Fa_N_top_m3 = np.swapaxes(Fa_N_top_swap, 0, 1)
    Fa_N_down_m3 = np.swapaxes(Fa_N_down_swap, 0, 1)

    return Fa_E_top_m3, Fa_E_down_m3, Fa_N_top_m3, Fa_N_down_m3


def get_stablefluxes(Fa_E, Fa_N, W):
    """Stabilize the outfluxes / influxes.

    During the reduced timestep the water cannot move further than 1/x * the
    gridcell, In other words at least x * the reduced timestep is needed to
    cross a gridcell.
    """
    Fa_E_abs = np.abs(Fa_E)
    Fa_N_abs = np.abs(Fa_N)

    stab = 1.0 / 2.0

    Fa_E_corrected = (Fa_E_abs / (Fa_E_abs + Fa_N_abs)) * stab * W[:-1, :, :]
    Fa_E_stable = np.minimum(Fa_E_abs, Fa_E_corrected)

    Fa_N_corrected = (Fa_N_abs / (Fa_E_abs + Fa_N_abs)) * stab * W[:-1, :, :]
    Fa_N_stable = np.minimum(Fa_N_abs, Fa_N_corrected)

    # get rid of the nan values
    Fa_E_stable[np.isnan(Fa_E_stable)] = 0
    Fa_N_stable[np.isnan(Fa_N_stable)] = 0

    # redefine
    Fa_E = np.sign(Fa_E) * Fa_E_stable
    Fa_N = np.sign(Fa_N) * Fa_N_stable

    return Fa_E, Fa_N


# Code
def getFa_Vert(
    Fa_E_top,
    Fa_E_down,
    Fa_N_top,
    Fa_N_down,
    E,
    P,
    W_top,
    W_down,
):

    # total moisture in the column
    W = W_top + W_down

    def divergence_zonal(Fa_E):
        """Define the horizontal fluxes over the boundaries."""

        # TODO: Is this correct? It looks like the eastern and western
        # boundaries are mixed up. Also the inserted zeros might cause trouble.
        # I think the implementation should be exactly like the meridional one.
        # I have verified that this code does exactly the same as the original.

        # fluxes over the eastern boundary
        Fa_E_boundary = np.zeros_like(Fa_E)
        Fa_E_boundary[:, :, :-1] = 0.5 * (Fa_E[:, :, :-1] + Fa_E[:, :, 1:])

        Fa_W_boundary = np.roll(Fa_E_boundary, 1)
        return Fa_W_boundary - Fa_E_boundary

    def divergence_meridional(Fa_N):
        """Define the horizontal fluxes over the boundaries."""
        Fa_N_boundary = np.zeros_like(Fa_N)
        Fa_N_boundary[:, 1:, :] = 0.5 * (Fa_N[:, :-1, :] + Fa_N[:, 1:, :])

        Fa_S_boundary = np.roll(Fa_N_boundary, -1, axis=1)
        return Fa_S_boundary - Fa_N_boundary

    zonal_divergence_top = divergence_zonal(Fa_E_top)
    zonal_divergence_down = divergence_zonal(Fa_E_down)
    meridional_divergence_top = divergence_meridional(Fa_N_top)
    meridional_divergence_down = divergence_meridional(Fa_N_down)

    # check the water balance
    residual_down = np.zeros_like(P)  # residual factor [m3]
    residual_top = np.zeros_like(P)  # residual factor [m3]

    tendency_down = (  # TODO why skip the N/S but not W/E edges?
        +zonal_divergence_down[:, 1:-1, :]
        + meridional_divergence_down[:, 1:-1, :]
        - P[:, 1:-1, :] * (W_down[:-1, 1:-1, :] / W[:-1, 1:-1, :])
        + E[:, 1:-1, :]
    )

    tendency_top = (  # TODO why skip the N/S but not W/E edges?
        +zonal_divergence_top[:, 1:-1, :]
        + meridional_divergence_top[:, 1:-1, :]
        - P[:, 1:-1, :] * (W_top[:-1, 1:-1, :] / W[:-1, 1:-1, :])
    )

    residual_down[:, 1:-1, :] = (
        W_down[1:, 1:-1, :] - W_down[:-1, 1:-1, :] - tendency_down
    )
    residual_top[:, 1:-1, :] = W_top[1:, 1:-1, :] - W_top[:-1, 1:-1, :] - tendency_top

    # compute the resulting vertical moisture flux; the vertical velocity so
    # that the new residual_down/W_down = residual_top/W_top (positive downward)
    Fa_Vert_raw = (
        W_down[1:, :, :] / W[1:, :, :] * (residual_down + residual_top) - residual_down
    )

    # stabilize the outfluxes / influxes; during the reduced timestep the
    # vertical flux can maximally empty/fill 1/x of the top or down storage
    stab = 1.0 / 4.0
    Fa_Vert_stable = np.minimum(
        np.abs(Fa_Vert_raw), np.minimum(stab * W_top[1:, :, :], stab * W_down[1:, :, :])
    )

    # redefine the vertical flux
    Fa_Vert = np.sign(Fa_Vert_raw) * Fa_Vert_stable

    return Fa_Vert
