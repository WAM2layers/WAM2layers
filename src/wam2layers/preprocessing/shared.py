"""Generic functions useful for preprocessing various input datasets."""

import numpy as np
import pandas as pd
import xarray as xr
from scipy.interpolate import interp1d


# old, keep only for reference and ecearth starting case
def resample(variable, divt, count_time, method="interp"):
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
        idx = int((t + 1) / divt + oddvector[t])

        if method == "bfill":
            new_var[t] = (1 / divt) * variable[idx - 1]
        elif method == "interp":
            new_var[t] = variable[idx - 1] + partvector[t] * (
                variable[idx] - variable[idx - 1]
            )
        else:
            raise ValueError(f"Unknown resample method {method}")

    return new_var


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

    a = (
        np.pi
        / 180.0
        * erad**2
        * grid_spacing
        * abs(np.sin(lat_s * np.pi / 180.0) - np.sin(lat_n * np.pi / 180.0))
    )

    # Calculate faces
    ly = grid_spacing * dg  # [m] length eastern/western boundary of a cell
    lx_n_gridcell = ly * np.cos((latitude + grid_spacing / 2) * np.pi / 180)
    lx_s_gridcell = ly * np.cos((latitude - grid_spacing / 2) * np.pi / 180)
    lx = 0.5 * (lx_n_gridcell + lx_s_gridcell)
    return a, ly, lx


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


def interpolate(x, xp, fp, axis=1, descending=False):
    """Linearly interpolate along an axis of an N-dimensional array.

    This function interpolates one slice at a time, i.e. if xp and fp are 4d
    arrays, x should be a 3d array and the function will return a 3d array.

    It is assumed that the input array is monotonic along the axis.
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
        ), "with desciending=True, xp must be monotonically increasing"

    # Check for out of bounds values
    if np.any(x[None, ...] < xp[0, ...]):
        raise ValueError("one or more x are below the lowest value of xp")
    if np.any(x[None, ...] > xp[-1, ...]):
        raise ValueError("one or more x are above the highest value of xp")

    # Find indices such that xp[lower] < x < xp[upper]
    upper = np.sum(x > xp, axis=0)
    lower = upper - 1

    # This will allow numpy advanced indexing to take an (N-1)D slice of an ND array
    upper = (upper, *np.meshgrid(*[range(l) for l in x.shape], indexing="ij"))
    lower = (lower, *np.meshgrid(*[range(l) for l in x.shape], indexing="ij"))

    fy = fp[lower] + (fp[upper] - fp[lower]) * (x - xp[lower]) / (xp[upper] - xp[lower])
    return fy


def insert_level(pressure_level_data, new_level, coord_value, dim_name="level"):
    """Insert a new level in the pressure level data.

    Note: new levels are inserted at the end of the dimension.
    Sorting by descending pressure is not performed automatically.

    Args:
        - pressure_level_data: xarray.DataArray with dimensions (time, lev, lat,
          lon)
        - new_level: the new data values that should be inserted as an
          additional level. If int or float, it will insert a constant value. If
          array like, it will insert the values in the array.
        - coord_value: coordinate value to use for the new level, e.g. 110000 Pa
          for well below the surface. Must be unique.
    """
    # Create dummy array that can be concatenated with the original data.
    dummy = pressure_level_data.isel({dim_name: 0}).copy()
    dummy[dim_name] = coord_value

    # Insert the nedim_name data into the dummy array
    if isinstance(new_level, xr.DataArray):
        dummy.values = new_level.values
    elif isinstance(new_level, np.ndarray):
        dummy.values = new_level
    elif isinstance(new_level, (int, float)):
        dummy.values = np.ones_like(dummy) * new_level
    else:
        raise ValueError("Invalid type for new_level")

    return xr.concat([pressure_level_data, dummy], dim=dim_name)


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


def change_units(data, target_freq):
    """Change units to m3.
    Multiply by edge length or area to get flux in m3
    Multiply by time to get accumulation instead of flux
    Divide by density of water to go from kg to m3
    """
    density = 1000  # [kg/m3]
    a, ly, lx = get_grid_info(data)

    total_seconds = pd.Timedelta(target_freq).total_seconds()

    for variable in data.data_vars:
        if variable in ["fx_upper", "fx_lower"]:
            data[variable] *= total_seconds / density * ly
        elif variable in ["fy_upper", "fy_lower"]:
            data[variable] *= total_seconds / density * lx[:, None]
        elif variable in ["evap", "precip"]:
            data[variable] *= total_seconds / density * a[:, None]
        elif variable in ["s_upper", "s_lower"]:
            data[variable] *= a[:, None] / density
        else:
            raise ValueError(f"Unrecognized variable {variable}")
        data[variable] = data[variable].assign_attrs(units="m**3")


def stabilize_fluxes(F, S, progress_tracker, config, t):
    """Stabilize the outfluxes / influxes.

    CFL: Water cannot move further than one grid cell per timestep.
    """
    coords = S.coords
    for level in ["upper", "lower"]:
        fx = F["fx_" + level]
        fy = F["fy_" + level]
        s = S["s_" + level]

        fx_abs = np.abs(fx)
        fy_abs = np.abs(fy)
        ft_abs = fx_abs + fy_abs

        # TODO: make 1/2 configurable?
        fx_limit = 1 / 2 * fx_abs / ft_abs * s.values
        fx_stable = np.minimum(fx_abs, fx_limit)

        fy_limit = 1 / 2 * fy_abs / ft_abs * s.values
        fy_stable = np.minimum(fy_abs, fy_limit)

        progress_tracker.track_stability_correction(fy_limit, fy_abs, config, coords, t)

        # Get rid of any nan values
        fx_stable.fillna(0)
        fy_stable.fillna(0)

        # Re-instate the sign
        F["fx_" + level] = np.sign(fx) * fx_stable
        F["fy_" + level] = np.sign(fy) * fy_stable


def convergence(fx, fy):
    # Note: latitude decreasing, hence positive fy gradient is convergence
    return np.gradient(fy, axis=-2) - np.gradient(fx, axis=-1)


def calculate_fz(F, S0, S1, kvf=0):
    """Calculate the vertical fluxes.

    The vertical flux is calculated as a closure term. Residuals are distributed
    proportionally to the amount of moisture in each layer.

    The flux is constrained such that it can never exceed

    Arguments:
        F: xarray dataset with fluxes evaluated at temporal midpoints between states
        S0: xarray dataset with states at current time t
        S1: xarray dataset with states at updated time t+1 (always forward looking)

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
        >>> calculate_fz(F, S0, S1)
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
    foo_upper = (
        (S1 - S0).s_upper
        - convergence(F.fx_upper, F.fy_upper)
        + F.precip.values * s_rel.s_upper
    )
    foo_lower = (
        (S1 - S0).s_lower
        - convergence(F.fx_lower, F.fy_lower)
        + F.precip.values * s_rel.s_lower
        - F.evap
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

    # Reinstate the sign
    return np.sign(fz) * fz_stable
