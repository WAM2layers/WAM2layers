from pathlib import Path

import click
import numpy as np
import pandas as pd
import xarray as xr
import yaml
from wam2layers.analysis.checks import check_input
from wam2layers.preprocessing.shared import (accumulation_to_flux,
                                                    calculate_humidity,
                                                    insert_level, interpolate,
                                                    sortby_ndarray)


def load_data(variable, date, config):
    """Load data for given variable and date."""
    prefix = config["filename_prefix"]
    # Include midnight of the next day (if available)
    extra = date + pd.Timedelta(days=1)

    if variable in ["u", "v", "q", "t"] and config["level_type"] == "model_levels":
        prefix += "_ml"

    filepath = Path(config["input_folder"]) / f"{prefix}_{variable}.nc"
    da = xr.open_dataset(filepath)[variable].sel(time=slice(date, extra))

    if "lev" in da.coords:
        da = da.rename(lev="level")

    if variable in ["u", "v", "q", "t"] and isinstance(config["levels"], list):
        return da.sel(level=config["levels"])

    return da.sel(time=slice(date, extra))


def preprocess_precip_and_evap(date, config):
    """Load and pre-process precipitation and evaporation."""
    # All incoming units are accumulations (in m) since previous time step
    evap = load_data("e", date, config)
    cp = load_data("cp", date, config)  # convective precipitation
    lsp = load_data("lsp", date, config)  # large scale precipitation
    precip = (cp + lsp)

    # Transfer negative (originally positive) values of evap to precip
    precip = np.maximum(precip, 0) + np.maximum(evap, 0)
    precip = precip

    # Change sign convention to all positive,
    evap = np.abs(np.minimum(evap, 0))

    precip = accumulation_to_flux(precip)
    evap = accumulation_to_flux(evap)

    return precip, evap


def get_edges(era5_modellevels):
    """Get the values of a and b at the edges of a subset of ERA5 modellevels."""
    # Load a and b coefficients
    table = Path(__file__).parent / "tableERA5model_to_pressure.csv"
    df = pd.read_csv(table)
    a = df["a [Pa]"].to_xarray().rename(index="level")
    b = df["b"].to_xarray().rename(index="level")

    # Calculate a and b at mid levels (model levels)
    a_full = ((a[1:] + a[:-1].values) / 2.0).sel(level=era5_modellevels)
    b_full = ((b[1:] + b[:-1].values) / 2.0).sel(level=era5_modellevels)

    # Interpolate to get parameters at edges between selected levels
    a_edge = xr.concat([a[0], (a_full[1:].values + a_full[:-1]) / 2, a[-1]], dim="level")
    b_edge = xr.concat([b[0], (b_full[1:].values + b_full[:-1]) / 2, b[-1]], dim="level")

    return a_edge, b_edge


def get_dp_modellevels(sp, modellevels):
    """Calculate pressure jump over subset of era5 modellevels."""
    if modellevels == "all":
        modellevels = list(range(138))

    sp_broadcast = sp.expand_dims({"level": modellevels}, axis=1)
    a_edge, b_edge = get_edges(modellevels)
    p_edges = a_edge + b_edge * sp_broadcast  # in Pa

    # calculate the difference between the pressure levels
    dp_modellevels = p_edges.diff(dim="level")  # in Pa

    return dp_modellevels


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
    for date in config["datelist"]:
        print(date)
        precip, evap = preprocess_precip_and_evap(date, config)

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

        # Integrate fluxes and states to upper and lower layer
        if config["level_type"] == "model_levels":
            # TODO: Check if this is a reasonable choice for boundary
            boundary = 111
            lower_layer = dp.level > boundary
            upper_layer = ~lower_layer
        if config["level_type"] == "pressure_levels":
            upper_layer = p < pb[:, None, :, :]
            lower_layer = pb[:, None, :, :] < p

        # Vertically integrate state over two layers
        s_lower = cw.where(lower_layer).sum(dim="level")
        s_upper = cw.where(upper_layer).sum(dim="level")

        # Determine the fluxes
        fx = u * cw  # eastward atmospheric moisture flux (kg m-1 s-1)
        fy = v * cw  # northward atmospheric moisture flux (kg m-1 s-1)

        # Vertically integrate fluxes over two layers
        fx_lower = fx.where(lower_layer).sum(dim="level")  # kg m-1 s-1
        fy_lower = fy.where(lower_layer).sum(dim="level")  # kg m-1 s-1
        fx_upper = fx.where(upper_layer).sum(dim="level")  # kg m-1 s-1
        fy_upper = fy.where(upper_layer).sum(dim="level")  # kg m-1 s-1

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
        check_input(ds)

        # Save preprocessed data
        filename = f"{date.strftime('%Y-%m-%d')}_fluxes_storages.nc"
        output_path = config["output_dir"] / filename
        ds.to_netcdf(output_path)

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
@click.argument('config_file', type=click.Path(exists=True))
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
