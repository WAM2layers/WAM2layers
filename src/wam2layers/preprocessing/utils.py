"""Generic functions useful for preprocessing various input datasets."""

from functools import lru_cache

import numpy as np
import pandas as pd
import xarray as xr
from scipy.interpolate import interp1d


def join_levels(pressure_level_data, surface_level_data):
    """Combine 3d pressure level and 2d surface level data.

    A dummy value of 1100 hPa is inserted for the surface pressure
    """
    bottom = pressure_level_data.isel(lev=0).copy()
    bottom.values = surface_level_data.values
    bottom["lev"] = 110000.0  # dummy value

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


def interpolate_old(old_var, old_pressure_levels, new_pressure_levels, type="linear"):
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


def interpolate(x, xp, fp, axis, descending=False) -> np.ndarray:
    """Linearly interpolate along an axis of an N-dimensional array.

    This function interpolates one slice at a time, i.e. if xp and fp are 4d
    arrays, x should be a 3d array and the function will return a 3d array.

    It is assumed that the input array is monotonic along the axis.

    Args:
        x: The coordinate to which you want to interpolate your data to.
        xp: The coordinates of the original data.
        fp: The original data.
        axis: The axis which should be interpolated over.
        descending: If the coordinates are in descending order (defaults to False)
    """
    # Cast input to numpy arrays
    x = np.asarray(x)
    xp = np.asarray(xp)
    fp = np.asarray(fp)

    # Move interpolation axis to first position for easier indexing
    xp = np.moveaxis(xp, axis, 0)
    fp = np.moveaxis(fp, axis, 0)

    # Handle descending axis
    if descending:
        xp = np.flip(xp, axis=0)
        fp = np.flip(fp, axis=0)
        assert (
            np.diff(xp, axis=0).min() >= 0
        ), "with descending=False, xp must be monotonically decreasing"
    else:
        assert (
            np.diff(xp, axis=0).min() >= 0
        ), "with descending=True, xp must be monotonically increasing"

    # Check for out of bounds values
    if np.any(x < xp[0, ...]):
        raise ValueError("one or more x are below the lowest value of xp")
    if np.any(x > xp[-1, ...]):
        raise ValueError("one or more x are above the highest value of xp")

    # Find indices such that xp[lower] < x < xp[upper]
    upper = np.sum(x > xp, axis=0)
    lower = upper - 1

    # This will allow numpy advanced indexing to take an (N-1)D slice of an ND array
    upper = (upper, *np.meshgrid(*[range(l) for l in x.shape], indexing="ij"))
    lower = (lower, *np.meshgrid(*[range(l) for l in x.shape], indexing="ij"))

    fy = fp[lower] + (fp[upper] - fp[lower]) * (x - xp[lower]) / (xp[upper] - xp[lower])
    return fy


def sortby_ndarray(array, other, axis):
    """Sort array along axis by the values in another array."""
    idx = np.argsort(other, axis=axis)
    return np.take_along_axis(array, idx, axis=axis)


def calculate_humidity(dewpoint, pressure):
    """
    Calculate the specific humidity from (surface) pressure and
    dew point temperature

    See further details at eq. 7.4 and 7.5 (Page 102) of:
    https://www.ecmwf.int/en/elibrary/20198-ifs-documentation-cy47r3-part-iv-physical-processes
    """
    Rd = 287.0597
    Rv = 461.5250
    a1 = 611.21
    a3 = 17.502
    a4 = 32.19
    t0 = 273.15

    # Calculation of saturation water vapour pressure from Teten's formula
    svp = a1 * np.exp(a3 * (dewpoint - t0) / (dewpoint - a4))

    # Specific humidity
    spec_hum = (Rd / Rv) * svp / (pressure - ((1 - Rd / Rv) * svp))

    return spec_hum


def accumulation_to_flux(data, input_frequency):
    """Convert precip and evap from accumulations to fluxes.

    Incoming data should have units of "m"

    Note: this does not currently modify the time labels. It can be argued
    that the times should be shifted. A good fix for this is welcome.
    """
    density = 1000  # [kg/m3]
    timestep = pd.Timedelta(input_frequency)
    nseconds = timestep.total_seconds()

    fluxdata = (density * data / nseconds).assign_attrs(units="kg m-2 s-1")

    # TODO: Adjust time points?
    # original_time = data.time
    # midpoint_time = original_time - timestep / 2
    # data["time"] = midpoint_time

    # Extrapolation introduces a small inconsistency at the last midnight...
    # data.interp(time=original_time, kwargs={"fill_value": "extrapolate"})
    return fluxdata


def add_bounds(ds: xr.Dataset) -> None:
    """Infer the lat and lon bounds and add to the dataset (in-place)."""
    lats = ds["latitude"].to_numpy()
    lons = ds["longitude"].to_numpy()
    res_lat = np.median((np.diff(lats)))  # infer resolution
    res_lon = np.median((np.diff(lons)))

    ds["latitude_bnds"] = (
        ("latitude", "bnds"),
        np.array((lats - res_lat / 2, lats + res_lat / 2)).T,
    )
    ds["longitude_bnds"] = (
        ("longitude", "bnds"),
        np.array((lons - res_lon / 2, lons + res_lon / 2)).T,
    )


def midpoints(x):
    """Linearly interpolate between the values of an array."""
    return (x[1:] + x[:-1]) / 2


@lru_cache(10)
def log_once(logger, msg: str):
    """Keep track of 10 different messages and then warn again.

    Adapted from: https://stackoverflow.com/a/66062313"""
    logger.info(msg)
