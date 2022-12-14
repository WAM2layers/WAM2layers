from functools import lru_cache
from pathlib import Path

import click
import numpy as np
import pandas as pd
import xarray as xr
import yaml

from wam2layers.preprocessing.shared import get_grid_info
from wam2layers.utils.profiling import Profiler, ProgressTracker


def parse_config(config_file):
    """Load and parse config file into dictionary."""
    with open(config_file) as f:
        config = yaml.safe_load(f)

    # Add datelist
    input_dates = pd.date_range(
        start=config["track_start_date"],
        end=config["track_end_date"],
        freq=config["input_frequency"],
        inclusive="left",
    )
    config["datelist"] = pd.date_range(
        start=input_dates[0],
        end=input_dates[-1],
        freq=config["target_frequency"],
        inclusive="right",
    )

    # Parse directories with pathlib
    config["preprocessed_data_folder"] = Path(
        config["preprocessed_data_folder"]
    ).expanduser()
    config["output_folder"] = Path(config["output_folder"]).expanduser() / "backtrack"

    # Check if input dir exists
    if not config["preprocessed_data_folder"].exists():
        raise ValueError(
            "Please create the preprocessed_data_folder before running the script"
        )

    # Create output dir if it doesn't exist yet
    if not config["output_folder"].exists():
        config["output_folder"].mkdir(parents=True)

    return config


def input_path(date, config):
    input_dir = config["preprocessed_data_folder"]
    return f"{input_dir}/{date.strftime('%Y-%m-%d')}_fluxes_storages.nc"


def output_path(date, config):
    output_dir = config["output_folder"]
    return f"{output_dir}/{date.strftime('%Y-%m-%d')}_s_track.nc"


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

    t1 = t.ceil(config["input_frequency"])
    da1 = read_data_at_time(t1)[variables[subset]]
    if t == t1:
        # this saves a lot of work if times are aligned with input
        return da1

    t0 = t.floor(config["input_frequency"])
    da0 = read_data_at_time(t0)[variables[subset]]
    if t == t0:
        return da0

    return da0 + (t - t0) / (t1 - t0) * (da1 - da0)


def load_fluxes(t):
    t_current = t - pd.Timedelta(config["target_frequency"]) / 2
    return load_data(t_current, "fluxes")


def load_states(t):
    t_prev = t
    t_next = t - pd.Timedelta(config["target_frequency"])
    states_prev = load_data(t_prev, "states")
    states_next = load_data(t_next, "states")
    return states_prev, states_next


def time_in_range(start, end, current):
    """Returns whether current is in the range [start, end]"""
    return start <= current <= end


def load_region(config):
    # TODO: make variable name more generic
    return xr.open_dataset(config["region"]).source_region


def to_edges_zonal(fx, periodic_boundary=False):
    """Define the horizontal fluxes over the east/west boundaries."""
    fe = np.zeros_like(fx)
    fe[:, :-1] = 0.5 * (fx[:, :-1] + fx[:, 1:])
    if periodic_boundary:
        fe[:, -1] = 0.5 * (fx[:, -1] + fx[:, 0])

    # find out where the positive and negative fluxes are
    f_pos = np.ones_like(fx)
    f_pos[fe < 0] = 0
    f_neg = f_pos - 1

    # separate directions west-east (all positive numbers)
    fe_we = fe * f_pos
    fe_ew = fe * f_neg

    # fluxes over the western boundary
    fw_we = look_west(fe_we)
    fw_ew = look_west(fe_ew)

    return fe_we, fe_ew, fw_we, fw_ew


def to_edges_meridional(fy):
    """Define the horizontal fluxes over the north/south boundaries."""
    fn = np.zeros_like(fy)
    fn[1:, :] = 0.5 * (fy[:-1, :] + fy[1:, :])

    # find out where the positive and negative fluxes are
    fn_pos = np.ones_like(fn)
    fn_pos[fn < 0] = 0  # invalid value encountered in less
    fn_neg = fn_pos - 1

    # separate directions south-north (all positive numbers)
    fn_sn = fn * fn_pos
    fn_ns = fn * fn_neg

    # fluxes over the southern boundary
    fs_sn = look_south(fn_sn)
    fs_ns = look_south(fn_ns)

    return fn_sn, fn_ns, fs_sn, fs_ns


def look_north(array):
    # Note: edges are reinserted at other end; but they're not used anyway
    return np.roll(array, 1, axis=-2)


def look_south(array):
    # Note: edges are reinserted at other end; but they're not used anyway
    return np.roll(array, -1, axis=-2)


def look_east(array):
    # Note: edges are reinserted at other end; but they're not used anyway
    return np.roll(array, -1, axis=-1)


def look_west(array):
    # Note: edges are reinserted at other end; but they're not used anyway
    return np.roll(array, 1, axis=-1)


def split_vertical_flux(Kvf, fv):
    f_downward = np.zeros_like(fv)
    f_upward = np.zeros_like(fv)
    f_downward[fv >= 0] = fv[fv >= 0]
    f_upward[fv <= 0] = fv[fv <= 0]
    f_upward = np.abs(f_upward)

    # include the vertical dispersion
    if Kvf != 0:
        f_upward = (1.0 + Kvf) * f_upward
        f_upward[fv >= 0] = fv[fv >= 0] * Kvf
        f_downward = (1.0 + Kvf) * f_downward
        f_downward[fv <= 0] = np.abs(fv[fv <= 0]) * Kvf

    return f_downward, f_upward


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


def stabilize_fluxes(current, previous):
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

        # Get rid of any nan values
        fx_stable.fillna(0)
        fy_stable.fillna(0)

        # Re-instate the sign
        current["fx_" + level] = np.sign(fx) * fx_stable
        current["fy_" + level] = np.sign(fy) * fy_stable


def convergence(fx, fy):
    # Note: latitude decreasing, hence positive fy gradient is convergence
    return np.gradient(fy, axis=-2) - np.gradient(fx, axis=-1)


def calculate_fv(fluxes, states_prev, states_next):
    """Calculate the vertical fluxes.

    Note: fluxes are given at temporal midpoints between states.
    """
    s_diff = states_prev - states_next
    s_mean = (states_prev + states_next) / 2
    s_total = s_mean.s_upper + s_mean.s_lower
    s_rel = s_mean / s_total

    tendency_upper = (
        convergence(fluxes.fx_upper, fluxes.fy_upper)
        - fluxes.precip.values * s_rel.s_upper
    )
    tendency_lower = (
        convergence(fluxes.fx_lower, fluxes.fy_lower)
        - fluxes.precip.values * s_rel.s_lower
        + fluxes.evap
    )

    residual_upper = s_diff.s_upper - tendency_upper
    residual_lower = s_diff.s_lower - tendency_lower

    # compute the resulting vertical moisture flux; the vertical velocity so
    # that the new residual_lower/s_lower = residual_upper/s_upper (positive downward)
    fv = s_rel.s_lower * (residual_upper + residual_lower) - residual_lower

    # stabilize the outfluxes / influxes; during the reduced timestep the
    # vertical flux can maximally empty/fill 1/x of the top or down storage
    stab = 1.0 / (config["kvf"] + 1.0)
    flux_limit = np.minimum(s_mean.s_upper, s_mean.s_lower)
    fv_stable = np.minimum(np.abs(fv), stab * flux_limit)

    # Reinstate the sign
    return np.sign(fv) * fv_stable


def backtrack(
    fluxes,
    states_prev,
    states_next,
    region,
    output,
):

    # Unpack input data
    fx_upper = fluxes["fx_upper"].values
    fy_upper = fluxes["fy_upper"].values
    fx_lower = fluxes["fx_lower"].values
    fy_lower = fluxes["fy_lower"].values
    evap = fluxes["evap"].values
    precip = fluxes["precip"].values
    f_vert = fluxes["f_vert"].values
    s_upper = states_prev["s_upper"].values
    s_lower = states_prev["s_lower"].values

    s_track_upper = output["s_track_upper_restart"].values
    s_track_lower = output["s_track_lower_restart"].values

    tagged_precip = region * precip
    s_total = s_upper + s_lower

    # separate the direction of the vertical flux and make it absolute
    f_downward, f_upward = split_vertical_flux(config["kvf"], f_vert)

    # Determine horizontal fluxes over the grid-cell boundaries
    f_e_lower_we, f_e_lower_ew, f_w_lower_we, f_w_lower_ew = to_edges_zonal(
        fx_lower, config["periodic_boundary"]
    )
    f_e_upper_we, f_e_upper_ew, f_w_upper_we, f_w_upper_ew = to_edges_zonal(
        fx_upper, config["periodic_boundary"]
    )

    (
        fy_n_lower_sn,
        fy_n_lower_ns,
        fy_s_lower_sn,
        fy_s_lower_ns,
    ) = to_edges_meridional(fy_lower)
    (
        fy_n_upper_sn,
        fy_n_upper_ns,
        fy_s_upper_sn,
        fy_s_upper_ns,
    ) = to_edges_meridional(fy_upper)

    # Short name for often used expressions
    s_track_relative_lower = s_track_lower / s_lower
    s_track_relative_upper = s_track_upper / s_upper
    if config["periodic_boundary"]:
        inner = np.s_[1:-1, :]
    else:
        inner = np.s_[1:-1, 1:-1]

    # Actual tracking (note: backtracking, all terms have been negated)
    s_track_lower[inner] += (
        +f_e_lower_we * look_east(s_track_relative_lower)
        + f_w_lower_ew * look_west(s_track_relative_lower)
        + fy_n_lower_sn * look_north(s_track_relative_lower)
        + fy_s_lower_ns * look_south(s_track_relative_lower)
        + f_upward * s_track_relative_upper
        - f_downward * s_track_relative_lower
        - fy_s_lower_sn * s_track_relative_lower
        - fy_n_lower_ns * s_track_relative_lower
        - f_e_lower_ew * s_track_relative_lower
        - f_w_lower_we * s_track_relative_lower
        + tagged_precip * (s_lower / s_total)
        - evap * s_track_relative_lower
    )[inner]

    s_track_upper[inner] += (
        +f_e_upper_we * look_east(s_track_relative_upper)
        + f_w_upper_ew * look_west(s_track_relative_upper)
        + fy_n_upper_sn * look_north(s_track_relative_upper)
        + fy_s_upper_ns * look_south(s_track_relative_upper)
        + f_downward * s_track_relative_lower
        - f_upward * s_track_relative_upper
        - fy_s_upper_sn * s_track_relative_upper
        - fy_n_upper_ns * s_track_relative_upper
        - f_w_upper_we * s_track_relative_upper
        - f_e_upper_ew * s_track_relative_upper
        + tagged_precip * (s_upper / s_total)
    )[inner]

    # down and top: redistribute unaccounted water that is otherwise lost from the sytem
    lower_to_upper = np.maximum(0, s_track_lower - states_next["s_lower"])
    upper_to_lower = np.maximum(0, s_track_upper - states_next["s_upper"])
    s_track_lower[inner] = (s_track_lower - lower_to_upper + upper_to_lower)[inner]
    s_track_upper[inner] = (s_track_upper - upper_to_lower + lower_to_upper)[inner]

    # Update output fields
    output["e_track"] += evap * (s_track_lower / s_lower)
    output["north_loss"] += (
        fy_n_upper_ns * s_track_relative_upper + fy_n_lower_ns * s_track_relative_lower
    )[1, :]
    output["south_loss"] += (
        fy_s_upper_sn * s_track_relative_upper + fy_s_lower_sn * s_track_relative_lower
    )[-2, :]

    if config["periodic_boundary"] == False:
        output["east_loss"] += (
            f_e_upper_ew * s_track_relative_upper
            + f_e_lower_ew * s_track_relative_lower
        )[:, -2]
        output["west_loss"] += (
            f_w_upper_we * s_track_relative_upper
            + f_w_lower_we * s_track_relative_lower
        )[:, 1]

    output["s_track_lower_restart"].values = s_track_lower
    output["s_track_upper_restart"].values = s_track_upper
    output["tagged_precip"] += tagged_precip


def initialize(config_file):
    """Read config, region, and initial states."""
    print(f"Initializing experiment with config file {config_file}")

    config = parse_config(config_file)
    region = load_region(config)

    output = initialize_outputs(region)
    region = region.values

    if config["restart"]:
        # Reload last state from existing output
        date = config["datelist"][-1] + pd.Timedelta(days=1)
        ds = xr.open_dataset(output_path(date, config))
        output["s_track_upper_restart"].values = ds.s_track_upper_restart.values
        output["s_track_lower_restart"].values = ds.s_track_lower_restart.values

    print(f"Output will be written to {config['output_folder']}.")
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
            "north_loss": xr.zeros_like(proto.isel(latitude=0, drop=True)),
            "south_loss": xr.zeros_like(proto.isel(latitude=0, drop=True)),
            "east_loss": xr.zeros_like(proto.isel(longitude=0, drop=True)),
            "west_loss": xr.zeros_like(proto.isel(longitude=0, drop=True)),
        }
    )
    return output


def write_output(output, t):
    # TODO: add back (and cleanup) coordinates and units
    path = output_path(t, config)
    print(f"{t} - Writing output to file {path}")
    output.astype("float32").to_netcdf(path)

    # Flush previous outputs
    output[
        [
            "e_track",
            "tagged_precip",
            "north_loss",
            "south_loss",
            "east_loss",
            "west_loss",
        ]
    ] *= 0


#############################################################################
# With the correct import statements, the code in the function below could
# alternatively be be used as a script in a separate python file or notebook.
#############################################################################


def run_experiment(config_file):
    """Run a backtracking experiment from start to finish."""
    global config

    config, region, output = initialize(config_file)

    progress_tracker = ProgressTracker(output)
    for t in reversed(config["datelist"]):

        fluxes = load_fluxes(t)
        states_prev, states_next = load_states(t)

        # Convert data to volumes
        change_units(states_prev, config["target_frequency"])
        change_units(states_next, config["target_frequency"])
        change_units(fluxes, config["target_frequency"])

        # Apply a stability correction if needed
        stabilize_fluxes(fluxes, states_next)

        # Determine the vertical moisture flux
        fluxes["f_vert"] = calculate_fv(fluxes, states_prev, states_next)

        # Only track the precipitation at certain dates
        if (
            time_in_range(
                config["event_start_date"],
                config["event_end_date"],
                t.strftime("%Y%m%d"),
            )
            == False
        ):
            fluxes["precip"] = 0

        backtrack(
            fluxes,
            states_prev,
            states_next,
            region,
            output,
        )

        # Daily output
        if t == t.floor(config["output_frequency"]):
            progress_tracker.print_progress(t, output)
            progress_tracker.store_intermediate_states(output)
            write_output(output, t)


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
    print("Welcome to WAM2layers.")
    print("Starting backtrack experiment.")
    run_experiment(config_file)


if __name__ == "__main__":
    cli()
