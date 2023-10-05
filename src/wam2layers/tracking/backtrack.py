import logging
from functools import lru_cache, partial

import click
import numpy as np
import pandas as pd
import xarray as xr

from wam2layers.config import Config
from wam2layers.preprocessing.shared import get_grid_info
from wam2layers.utils import load_region
from wam2layers.utils.profiling import ProgressTracker

logger = logging.getLogger(__name__)


def get_tracking_dates(config):
    """Dates for tracking."""
    # E.g. if data from 00:00 to 23:00
    # Using target freq directly would end at 23:45 or so
    input_dates = pd.date_range(
        start=config.track_start_date,
        end=config.track_end_date,
        freq=config.input_frequency,
    )

    return pd.date_range(
        start=input_dates[0],
        end=input_dates[-1],
        freq=config.target_frequency,
        inclusive="right",
    )


def input_path(date, config):
    input_dir = config.preprocessed_data_folder
    return f"{input_dir}/{date.strftime('%Y-%m-%d')}_fluxes_storages.nc"


def output_path(date, config):
    output_dir = config.output_folder
    return f"{output_dir}/backtrack_{date.strftime('%Y-%m-%dT%H-%M')}.nc"


# LRU Cache keeps the file open so we save a bit on I/O
@lru_cache(maxsize=2)
def read_data_at_date(d):
    """Load input data for given date."""
    file = input_path(d, config)
    return xr.open_dataset(file, cache=False)


# This one can't be cached as we'll be overwriting the content.
def read_data_at_time(t):
    """Get a single time slice from input data at time t."""
    ds = read_data_at_date(t)
    return ds.sel(time=t, drop=True)


def load_data(t, subset="fluxes"):
    """Load variable at t, interpolate if needed."""
    variables = {
        "fluxes": ["fx_upper", "fx_lower", "fy_upper", "fy_lower", "evap", "precip"],
        "states": ["s_upper", "s_lower"],
    }

    t1 = t.ceil(config.input_frequency)
    da1 = read_data_at_time(t1)[variables[subset]]
    if t == t1:
        # this saves a lot of work if times are aligned with input
        return da1

    t0 = t.floor(config.input_frequency)
    da0 = read_data_at_time(t0)[variables[subset]]
    if t == t0:
        return da0

    return da0 + (t - t0) / (t1 - t0) * (da1 - da0)


def load_region(config):
    # TODO: make variable name more generic
    return xr.open_dataset(config.region).source_region


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

    vq = np.where(vp > 0, vp * qp[south, inner], vp * qp[north, inner])  # [M, N+2]

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


def vertical_dispersion(fv, s_lower, s_upper):
    """Calculate additional vertical mixing due to convective dispersion.

    dispersion = kvf * |Fv| * dS/dz
    """
    return config.kvf * np.abs(fv) * (s_upper - s_lower)


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


def stabilize_fluxes(current, previous, progress_tracker, t):
    """Stabilize the outfluxes / influxes.

    CFL: Water cannot move further than one grid cell per timestep.
    """
    for level in ["upper", "lower"]:
        fx = current["fx_" + level]
        fy = current["fy_" + level]
        s = previous["s_" + level]

        fx_abs = np.abs(fx)
        fy_abs = np.abs(fy)
        ft_abs = fx_abs + fy_abs

        fx_corrected = 1 / 2 * fx_abs / ft_abs * s.values
        fx_stable = np.minimum(fx_abs, fx_corrected)

        fy_corrected = 1 / 2 * fy_abs / ft_abs * s.values
        fy_stable = np.minimum(fy_abs, fy_corrected)

        progress_tracker.track_stability_correction(fy_corrected, fy_abs, config, t)

        # Get rid of any nan values
        fx_stable.fillna(0)
        fy_stable.fillna(0)

        # Re-instate the sign
        current["fx_" + level] = np.sign(fx) * fx_stable
        current["fy_" + level] = np.sign(fy) * fy_stable


def convergence(fx, fy):
    # Note: latitude decreasing, hence positive fy gradient is convergence
    return np.gradient(fy, axis=-2) - np.gradient(fx, axis=-1)


def calculate_fz(F, S0, S1):
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
    """
    s_mean = (S1 + S0) / 2
    s_total = s_mean.s_upper + s_mean.s_lower
    s_rel = s_mean / s_total

    residual_upper = (
        (S1 - S0).s_upper
        - convergence(F.fx_upper, F.fy_upper)
        + F.precip.values * s_rel.s_upper
    )
    residual_lower = (
        (S1 - S0).s_lower
        - convergence(F.fx_lower, F.fy_lower)
        + F.precip.values * s_rel.s_lower
        - F.evap
    )

    # compute the resulting vertical moisture flux; the vertical velocity so
    # that the new residual_lower/s_lower = residual_upper/s_upper (positive downward)
    fz = s_rel.s_lower * (residual_upper + residual_lower) - residual_lower

    # stabilize the outfluxes / influxes; during the reduced timestep the
    # vertical flux can maximally empty/fill 1/x of the top or down storage
    stab = 1.0 / (config.kvf + 1.0)
    flux_limit = np.minimum(
        s_mean.s_upper, s_mean.s_lower
    )  # TODO why is this not 'just' the upstream bucket?
    fz_stable = np.minimum(np.abs(fz), stab * flux_limit)

    # Reinstate the sign
    return np.sign(fz) * fz_stable


def backtrack(
    F,
    S1,
    S0,
    region,
    output,
):
    # Unpack input data
    # TODO move staggering to preprocessing (?)
    fx_upper = stagger_x(F["fx_upper"].values)
    fy_upper = stagger_y(F["fy_upper"].values)
    fx_lower = stagger_x(F["fx_lower"].values)
    fy_lower = stagger_y(F["fy_lower"].values)
    evap = F["evap"].values
    precip = F["precip"].values
    f_vert = F["f_vert"].values
    s_upper = S1["s_upper"].values
    s_lower = S1["s_lower"].values

    s_track_upper = output["s_track_upper_restart"].values
    s_track_lower = output["s_track_lower_restart"].values

    tagged_precip = region * precip
    s_total = s_upper + s_lower

    # Determine fluxes at the faces of the grid cells

    # Short name for often used expressions
    s_track_relative_lower = np.minimum(s_track_lower / s_lower, 1.0)
    s_track_relative_upper = np.minimum(s_track_upper / s_upper, 1.0)

    # Actual tracking (note: backtracking, all fluxes change sign)
    bc = config.periodic_boundary
    s_track_lower += (
        +horizontal_advection(s_track_relative_lower, -fx_lower, -fy_lower, bc)
        + vertical_advection(-f_vert, s_track_relative_lower, s_track_relative_upper)
        + vertical_dispersion(-f_vert, s_track_relative_lower, s_track_relative_upper)
        + (tagged_precip * (s_lower / s_total))
        - (evap * s_track_relative_lower)
    )

    s_track_upper += (
        +horizontal_advection(s_track_relative_upper, -fx_upper, -fy_upper, bc)
        - vertical_advection(-f_vert, s_track_relative_lower, s_track_relative_upper)
        - vertical_dispersion(-f_vert, s_track_relative_lower, s_track_relative_upper)
        + (tagged_precip * (s_upper / s_total))
    )

    # down and top: redistribute unaccounted water that is otherwise lost from the sytem
    lower_to_upper = np.maximum(0, s_track_lower - S0["s_lower"])
    upper_to_lower = np.maximum(0, s_track_upper - S0["s_upper"])
    s_track_lower = s_track_lower - lower_to_upper + upper_to_lower
    s_track_upper = s_track_upper - upper_to_lower + lower_to_upper

    # Update output fields
    output["e_track"] += evap * np.minimum(s_track_lower / s_lower, 1.0)

    # Bookkeep boundary losses as "tracked moisture at grid edges"
    output["e_track"][0, :] += (s_track_upper + s_track_lower)[0, :]
    output["e_track"][-1, :] += (s_track_upper + s_track_lower)[-1, :]
    s_track_upper[0, :] = 0
    s_track_upper[-1, :] = 0
    s_track_lower[0, :] = 0
    s_track_lower[-1, :] = 0
    if config.periodic_boundary is False:
        output["e_track"][:, 0] += (s_track_upper + s_track_lower)[:, 0]
        output["e_track"][:, -1] += (s_track_upper + s_track_lower)[:, -1]
        s_track_upper[:, 0] = 0
        s_track_upper[:, -1] = 0
        s_track_lower[:, 0] = 0
        s_track_lower[:, -1] = 0

    output["s_track_lower_restart"].values = s_track_lower
    output["s_track_upper_restart"].values = s_track_upper
    output["tagged_precip"] += tagged_precip


def initialize(config_file):
    """Read config, region, and initial states."""
    logger.info(f"Initializing experiment with config file {config_file}")

    config = Config.from_yaml(config_file)
    region = load_region(config)

    output = initialize_outputs(region)
    region = region.values

    if config.restart:
        # Reload last state from existing output
        tracking_dates = get_tracking_dates(config)
        date = tracking_dates[-1] + pd.Timedelta(days=1)
        ds = xr.open_dataset(output_path(date, config))
        output["s_track_upper_restart"].values = ds.s_track_upper_restart.values
        output["s_track_lower_restart"].values = ds.s_track_lower_restart.values

    logger.info(f"Output will be written to {config.output_folder}.")
    return config, region, output


def initialize_outputs(region):
    """Allocate output arrays."""

    proto = region
    output = xr.Dataset(
        {
            # Keep last state for a restart
            "s_track_upper_restart": xr.zeros_like(proto),
            "s_track_lower_restart": xr.zeros_like(proto),
            "e_track": xr.zeros_like(proto),
            "tagged_precip": xr.zeros_like(proto),
        }
    )
    return output


def write_output(output, t):
    # TODO: add back (and cleanup) coordinates and units
    path = output_path(t, config)
    logger.info(f"{t} - Writing output to file {path}")
    output.astype("float32").to_netcdf(path)

    # Flush previous outputs
    output[["e_track", "tagged_precip"]] *= 0


#############################################################################
# With the correct import statements, the code in the function below could
# alternatively be be used as a script in a separate python file or notebook.
#############################################################################


def run_experiment(config_file):
    """Run a backtracking experiment from start to finish."""
    global config

    config, region, output = initialize(config_file)

    tracking_dates = get_tracking_dates(config)
    event_start, event_end = config.event_start_date, config.event_end_date

    progress_tracker = ProgressTracker(output)
    for t1 in reversed(tracking_dates):
        t0 = t1 - pd.Timedelta(config.target_frequency)
        th = t1 - pd.Timedelta(config.target_frequency) / 2

        S1 = load_data(t1, "states")
        S0 = load_data(t0, "states")
        F = load_data(th, "fluxes")

        # Convert data to volumes
        change_units(S1, config.target_frequency)
        change_units(S0, config.target_frequency)
        change_units(F, config.target_frequency)

        # Apply a stability correction if needed
        stabilize_fluxes(F, S0, progress_tracker, t1)

        # Determine the vertical moisture flux
        F["f_vert"] = calculate_fz(F, S0, S1)

        # Only track the precipitation at certain timesteps
        if not event_start <= t1 <= event_end:
            F["precip"] = 0

        backtrack(F, S1, S0, region, output)

        # Daily output
        if t1 == t1.floor(config.output_frequency) or t1 == tracking_dates[0]:
            progress_tracker.print_progress(t1, output)
            progress_tracker.store_intermediate_states(output)
            write_output(output, t1)

    logger.info("Experiment complete.")


###########################################################################
# The code below makes it possible to run wam2layers from the command line:
# >>> python backtrack.py path/to/cases/era5_2021.yaml
# or even:
# >>> wam2layers backtrack path/to/cases/era5_2021.yaml
###########################################################################


@click.command()
@click.argument("config_file", type=click.Path(exists=True))
def cli(config_file):
    """Run WAM2layers backtrack experiment from the command line.

    CONFIG_FILE: Path to WAM2layers experiment configuration file.

    Usage examples:

        \b
        - python path/to/backtrack.py path/to/cases/era5_2021.yaml
        - wam2layers backtrack path/to/cases/era5_2021.yaml
    """
    logger.info("Welcome to WAM2layers.")
    logger.info("Starting backtrack experiment.")
    run_experiment(config_file)


if __name__ == "__main__":
    cli()
