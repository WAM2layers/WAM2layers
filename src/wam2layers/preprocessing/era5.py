from pathlib import Path

import click
import numpy as np
import pandas as pd
import xarray as xr
import yaml

from wam2layers.analysis.checks import check_input
from wam2layers.preprocessing.shared import (accumulation_to_flux,
                                             calculate_humidity, insert_level,
                                             interpolate, sortby_ndarray)


class InvalidConfig(Exception):
    pass


def load_data(variable, date, config):
    """Load data for given variable and date."""
    template = config["filename_template"]

    # If it's a 4d variable, we need to set the level type
    if variable in ["u", "v", "q"]:
        if config["level_type"] == "model_levels":
            prefix = "_ml"
        elif config["level_type"] == "pressure_levels":
            prefix = "_pl"
    # If it's a 3d variable we do not need to set the level type
    else:
        prefix = ""

    # Load data
    filepath = template.format(
        year=date.year,
        month=date.month,
        day=date.day,
        levtype=prefix,
        variable=variable,
    )
    da = xr.open_dataset(filepath, chunks=config["chunks"]).sel(
        time=date.strftime("%Y%m%d")
    )[variable]

    if "lev" in da.coords:
        da = da.rename(lev="level")

    # If it's 4d data we want to select a subset of the levels
    if variable in ["u", "v", "q"] and isinstance(config["levels"], list):
        return da.sel(level=config["levels"])

    return da


def preprocess_precip_and_evap(date, config):
    """Load and pre-process precipitation and evaporation."""
    # All incoming units are accumulations (in m) since previous time step
    evap = load_data("e", date, config)
    precip = load_data("tp", date, config)  # total precipitation

    # Transfer negative (originally positive) values of evap to precip
    precip = np.maximum(precip, 0) + np.maximum(evap, 0)
    precip = precip

    # Change sign convention to all positive,
    evap = np.abs(np.minimum(evap, 0))

    precip = accumulation_to_flux(precip)
    evap = accumulation_to_flux(evap)
    return precip, evap


def load_era5_ab():
    """Load a and b coefficients."""
    table = Path(__file__).parent / "tableERA5model_to_pressure.csv"
    df = pd.read_csv(table)
    a = df["a [Pa]"].values
    b = df.b.values
    return a, b


def midpoints(x):
    """Linearly interpolate between the values of an array."""
    return (x[1:] + x[:-1]) / 2


def get_edges(era5_modellevels):
    """Get the values of a and b at the edges of a subset of ERA5 modellevels."""
    a, b = load_era5_ab()

    if era5_modellevels == "all":
        return a, b

    # With python indexing starting at 0 and 137 full levels,
    # count from 0-136 instead of 1-137.
    era5_modellevels = [i - 1 for i in era5_modellevels]

    # Calculate the values of a and be at midpoints and extract levels
    a_full = midpoints(a)[era5_modellevels]
    b_full = midpoints(b)[era5_modellevels]

    # Calculate the values of a and b at the edges *between* these levels
    a_edge = midpoints(a_full)
    b_edge = midpoints(b_full)

    # Reinsert the original outer edges
    a_edge = np.block([a[0], a_edge, a[-1]])
    b_edge = np.block([b[0], b_edge, b[-1]])

    return a_edge, b_edge


def get_dp_modellevels(sp, levels):
    """Calculate pressure jump over subset of era5 modellevels."""
    a, b = get_edges(levels)

    # Use dummy level -1 to avoid collision with existing levels
    # first coord will be discarded anyway on taking the diff below.
    a = xr.DataArray(a, coords={"level": np.insert(levels, 0, -1)})
    b = xr.DataArray(b, coords={"level": np.insert(levels, 0, -1)})

    p_edge = a + b * sp

    # Calculate the difference between the pressure levels
    dp = p_edge.diff("level")  # in Pa

    return dp


def get_dp_pressurelevels(q, u, v, ps, qs, us, vs):
    """Get dp with consistent u, v, q for ERA5 pressure level data."""
    p = u.level.broadcast_like(u) * 100  # Pa

    # Insert top of atmosphere values
    # Assume wind at top same as values at lowest pressure, humidity at top 0
    u = insert_level(u, u.isel(level=0), 0)
    v = insert_level(v, v.isel(level=0), 0)
    q = insert_level(q, 0, 0)
    p = insert_level(p, 0, 0)

    # Insert surface level values (at a high dummy pressure value)
    u = insert_level(u, us, 110000)
    v = insert_level(v, vs, 110000)
    q = insert_level(q, qs, 110000)
    p = insert_level(p, ps, 110000)

    # Sort arrays by pressure (ascending)
    u.values = sortby_ndarray(u.values, p.values, axis=1)
    v.values = sortby_ndarray(v.values, p.values, axis=1)
    q.values = sortby_ndarray(q.values, p.values, axis=1)
    p.values = sortby_ndarray(p.values, p.values, axis=1)

    # Insert boundary level values (at a ridiculous dummy pressure value)
    p_boundary = 0.72878581 * np.array(ps) + 7438.803223
    u = insert_level(u, interpolate(p_boundary, p, u), 150000)
    v = insert_level(v, interpolate(p_boundary, p, v), 150000)
    q = insert_level(q, interpolate(p_boundary, p, q), 150000)
    p = insert_level(p, p_boundary, 150000)

    # Sort arrays by pressure once more (ascending)
    u.values = sortby_ndarray(u.values, p.values, axis=1)
    v.values = sortby_ndarray(v.values, p.values, axis=1)
    q.values = sortby_ndarray(q.values, p.values, axis=1)
    p.values = sortby_ndarray(p.values, p.values, axis=1)

    # Reset level coordinate as its values have become meaningless
    nlev = u.level.size
    u = u.assign_coords(level=np.arange(nlev))
    v = v.assign_coords(level=np.arange(nlev))
    q = q.assign_coords(level=np.arange(nlev))
    p = p.assign_coords(level=np.arange(nlev))

    # Calculate pressure jump
    dp = p.diff("level")
    assert np.all(dp > 0), "Pressure levels should increase monotonically"

    # Interpolate to midpoints
    midpoints = 0.5 * (u.level.values[1:] + u.level.values[:-1])
    dp = dp.assign_coords(level=midpoints)
    u = u.interp(level=midpoints)
    v = v.interp(level=midpoints)
    q = q.interp(level=midpoints)
    p = p.interp(level=midpoints)

    # mask values below surface
    above_surface = p < np.array(ps)[:, None, :, :]
    u = u.where(above_surface)
    v = v.where(above_surface)
    q = q.where(above_surface)
    p = p.where(above_surface)

    return dp, p, q, u, v, p_boundary


def parse_config(config_file):
    """Read and parse case configuration file."""
    with open(config_file) as f:
        config = yaml.safe_load(f)

    # Create the preprocessed data folder if it does not exist yet
    config["output_dir"] = Path(config["preprocessed_data_folder"]).expanduser()
    config["output_dir"].mkdir(exist_ok=True, parents=True)

    config["datelist"] = pd.date_range(
        start=config["preprocess_start_date"],
        end=config["preprocess_end_date"],
        freq="d",
        inclusive="left",
    )
    return config


def prep_experiment(config_file):
    """Pre-process all data for a given config file.

    This function expects the following configuration settings:

    - preprocess_start_date: formatted as YYYYMMDD, e.g. '20210701'
    - preprocess_end_date: formatted as YYYYMMDD, e.g. '20210716'
    - level_type: either "pressure_levels" or "model_levels"
    - levels: "all" or a list of integers with the desired (model or pressure)
      levels.
    - input_folder: path where raw era5 input data can be found, e.g.
      /home/peter/WAM2layers/era5_2021
    - preprocessed_data_folder: path where preprocessed data should be stored.
      This directory will be created if it does not exist. E.g.
      /home/peter/WAM2layers/preprocessed_data_2021
    - filename_prefix: Fixed part of filename. This function will infer the
      variable name and add _ml for model level data. E.g. with prefix =
      "FloodCase_202107" this function will be able to find
      FloodCase_202107_ml_u.nc or FloodCase_202107_u.nc and
      FloodCase_202107_sp.nc
    """
    config = parse_config(config_file)

    if config["chunks"] is not None:
        print("Starting dask cluster")
        from dask.distributed import Client

        client = Client()
        print(f"To see the dask dashboard, go to {client.dashboard_link}")

    for date in config["datelist"]:
        print(date)

        # 4d fields
        levels = config["levels"]
        q = load_data("q", date, config)  # in kg kg-1
        u = load_data("u", date, config)  # in m/s
        v = load_data("v", date, config)  # in m/s
        sp = load_data("sp", date, config)  # in Pa

        if config["level_type"] == "model_levels":
            dp = get_dp_modellevels(sp, levels)

        if config["level_type"] == "pressure_levels":
            d2m = load_data("d2m", date, config)  # Dew point in K
            q2m = calculate_humidity(d2m, sp)  # kg kg-1
            u10 = load_data("u10", date, config)  # in m/s
            v10 = load_data("v10", date, config)  # in m/s
            dp, p, q, u, v, pb = get_dp_pressurelevels(q, u, v, sp, q2m, u10, v10)

        # Calculate column water vapour
        g = 9.80665  # gravitational accelleration [m/s2]
        cwv = q * dp / g  # (kg/m2)

        try:
            # Calculate column water instead of column water vapour
            tcw = load_data("tcw", date, config)  # kg/m2
            cw = (tcw / cwv.sum(dim="level")) * cwv  # column water (kg/m2)
            # TODO: warning if cw >> cwv
        except FileNotFoundError:
            # Fluxes will be calculated based on the column water vapour
            cw = cwv

        # Determine the fluxes
        fx = u * cw  # eastward atmospheric moisture flux (kg m-1 s-1)
        fy = v * cw  # northward atmospheric moisture flux (kg m-1 s-1)

        # Integrate fluxes and states to upper and lower layer
        if config["level_type"] == "model_levels":
            # TODO: Check if this is a reasonable choice for boundary
            boundary = 111
            idx = dp.level.searchsorted(boundary, side="right")
            upper = np.s_[:, :idx, :, :]
            lower = np.s_[:, idx:, :, :]
            s_lower = cw[lower].sum(dim="level")
            s_upper = cw[upper].sum(dim="level")
            fx_lower = fx[lower].sum(dim="level")  # kg m-1 s-1
            fy_lower = fy[lower].sum(dim="level")  # kg m-1 s-1
            fx_upper = fx[upper].sum(dim="level")  # kg m-1 s-1
            fy_upper = fy[upper].sum(dim="level")  # kg m-1 s-1

        elif config["level_type"] == "pressure_levels":
            upper_layer = p < pb[:, None, :, :]
            lower_layer = pb[:, None, :, :] < p

            # Vertically integrate state over two layers
            s_lower = cw.where(lower_layer).sum(dim="level")
            s_upper = cw.where(upper_layer).sum(dim="level")
            fx_lower = fx.where(lower_layer).sum(dim="level")  # kg m-1 s-1
            fy_lower = fy.where(lower_layer).sum(dim="level")  # kg m-1 s-1
            fx_upper = fx.where(upper_layer).sum(dim="level")  # kg m-1 s-1
            fy_upper = fy.where(upper_layer).sum(dim="level")  # kg m-1 s-1

        else:
            raise InvalidConfig(
                "Expected config level_type to be one of ['model_level', "
                f"'pressure_level'], but got {config['level_type']}."
            )

        precip, evap = preprocess_precip_and_evap(date, config)

        # Combine everything into one dataset
        ds = xr.Dataset(
            {
                "fx_upper": fx_upper.assign_attrs(units="kg m-1 s-1"),
                "fy_upper": fy_upper.assign_attrs(units="kg m-1 s-1"),
                "fx_lower": fx_lower.assign_attrs(units="kg m-1 s-1"),
                "fy_lower": fy_lower.assign_attrs(units="kg m-1 s-1"),
                "s_upper": s_upper.assign_attrs(units="kg m-2"),
                "s_lower": s_lower.assign_attrs(units="kg m-2"),
                "evap": evap,
                "precip": precip,
            }
        )

        # Verify that the data meets all the requirements for the model
        # check_input(ds)

        # Save preprocessed data
        filename = f"{date.strftime('%Y-%m-%d')}_fluxes_storages.nc"
        output_path = config["output_dir"] / filename
        ds.astype("float32").to_netcdf(output_path)

    # Close the dask cluster when done
    if config["chunks"] is not None:
        client.shutdown()


################################################################################
# To run this script interactively in e.g. Spyder, uncomment the following line:
# prep_experiment("../../cases/era5_2021.yaml")
################################################################################

###########################################################################
# The code below makes it possible to run wam2layers from the command line:
# >>> python backtrack.py path/to/cases/era5_2021.yaml
# or even:
# >>> wam2layers backtrack path/to/cases/era5_2021.yaml
###########################################################################


@click.command()
@click.argument("config_file", type=click.Path(exists=True))
def cli(config_file):
    """Preprocess ERA5 data for WAM2layers tracking experiments.

    CONFIG_FILE: Path to WAM2layers experiment configuration file.

    Usage examples:

        \b
        - python path/to/preprocessing/era5.py path/to/cases/era5_2021.yaml
        - wam2layers preprocess era5 path/to/cases/era5_2021.yaml
    """
    prep_experiment(config_file)


if __name__ == "__main__":
    cli()
