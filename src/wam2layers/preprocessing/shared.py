"""General preprocessing routine."""
import pandas as pd
import xarray as xr

from wam2layers.config import Config
from wam2layers.preprocessing import era5, logger
from wam2layers.preprocessing.input_validation import validate_input
from wam2layers.preprocessing.pressure_levels import (
    extend_pressurelevels,
    interp_dp_midpoints,
    mask_below_surface_data,
)
from wam2layers.preprocessing.utils import add_bounds
from wam2layers.preprocessing.xarray_append import append_to_netcdf
from wam2layers.reference.variables import PREPROCESSED_DATA_ATTRIBUTES


def get_input_dates(config: Config) -> pd.DatetimeIndex:
    """Dates for pre-processing."""
    return pd.date_range(
        start=config.preprocess_start_date,
        end=config.preprocess_end_date,
        freq=config.input_frequency,
    )


def get_input_data(
    datetime: pd.Timestamp, config: Config, source: str
) -> tuple[xr.Dataset, dict[str, str]]:
    """Retrieve input data from a specified source (e.g. ERA5, CMIP6).

    This will load the xr.Dataset required for the preprocessing routines.
        This dataset has the following coordinates: latitude, longitude, level.
            Note that for pressure level data, the level coordinates should be the air
            pressure (Pa). For model level data these are integers, counting down
            towards the surface.
        It has the following variables:
            q: Specific humidity at pressure levels (kg/kg)
            u: Eastward horizontal wind speed at pressure levels (m/s)
            v: Northward horizontal wind speed at pressure levels (m/s)
            ps: Air pressure at the surface (Pa)
            twc (optional): Total column water content (kg/m2)
        If the data uses model layers, this dataset has the variable:
            dp: Pressure difference over each model layer (Pa).
        If the data uses pressure layers, the dataset will have the following additional
        variables:
            qs: Specific humidity at the surface (kg/kg)
            us: Eastward horizontal wind speed at the surface (m/s)
            vs: Northward horizontal wind speed at the surface (m/s)

    Args:
        datetime: Which day of data should be loaded
        config: WAM2layers configuration
        source (optional): Which source to use. Defaults to "ERA5".

    Returns:
        Prepared input data
        Attributes corresponding to the input data source
    """
    if source == "ERA5":
        data, input_data_attrs = era5.get_input_data(datetime, config)
    else:
        msg = "Only the ERA5 input data loader has been implemented."
        raise NotImplementedError(msg)

    validate_input(data, config)
    return data, input_data_attrs


def prep_experiment(config_file: Config, data_source: str):
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
    config = Config.from_yaml(config_file)

    for datetime in get_input_dates(config):
        is_new_day = datetime == datetime.floor("1d")
        logger.info(datetime)

        input_data, input_data_attrs = get_input_data(datetime, config, data_source)
        surface_data = input_data.drop_dims("level")
        level_data = input_data.drop_vars(surface_data.data_vars)

        if config.level_type == "pressure_levels":
            level_data, pb = extend_pressurelevels(level_data, surface_data, config)
            level_data, dp = interp_dp_midpoints(level_data)
            level_data = mask_below_surface_data(level_data, surface_data["ps"])
        else:
            dp = level_data["dp"]

        # Calculate column water vapour
        g = 9.80665  # gravitational acceleration [m/s2]
        cwv = level_data["q"] * dp / g  # (kg/m2)

        if "tcw" in surface_data:
            # Calculate column water instead of column water vapour
            correction = surface_data["tcw"] / cwv.sum(dim="level")
            cw = correction * cwv  # column water (kg/m2)
            if is_new_day:
                logger.info(
                    "Total column water correction:\n"
                    "    ratio total column water / computed integrated water vapour\n"
                    f"    mean over grid for this timestep {correction.mean().item():.4f}"
                )
        else:  # Fluxes will be calculated based on the column water vapour
            cw = cwv

        # Integrate fluxes and states to upper and lower layer
        if config.level_type == "model_levels":
            lower_layer = dp["level"] > config.level_layer_boundary
            upper_layer = ~lower_layer

        if config.level_type == "pressure_levels":
            upper_layer = level_data["p"] < pb.broadcast_like(level_data["p"])
            lower_layer = ~upper_layer

        # Vertically integrate state over two layers
        s_lower = cw.where(lower_layer).sum(dim="level")
        s_upper = cw.where(upper_layer).sum(dim="level")

        # Determine the fluxes
        fx = level_data["u"] * cw  # eastward atmospheric moisture flux (kg m-1 s-1)
        fy = level_data["v"] * cw  # northward atmospheric moisture flux (kg m-1 s-1)

        # Vertically integrate fluxes over two layers
        fx_lower = fx.where(lower_layer).sum(dim="level")  # kg m-1 s-1
        fy_lower = fy.where(lower_layer).sum(dim="level")  # kg m-1 s-1
        fx_upper = fx.where(upper_layer).sum(dim="level")  # kg m-1 s-1
        fy_upper = fy.where(upper_layer).sum(dim="level")  # kg m-1 s-1

        # Combine everything into one dataset
        ds = (
            xr.Dataset(
                {
                    "fx_upper": fx_upper.assign_attrs(units="kg m-1 s-1"),
                    "fy_upper": fy_upper.assign_attrs(units="kg m-1 s-1"),
                    "fx_lower": fx_lower.assign_attrs(units="kg m-1 s-1"),
                    "fy_lower": fy_lower.assign_attrs(units="kg m-1 s-1"),
                    "s_upper": s_upper.assign_attrs(units="kg m-2"),
                    "s_lower": s_lower.assign_attrs(units="kg m-2"),
                    "evap": surface_data["evap"],
                    "precip": surface_data["precip"],
                }
            )
            .expand_dims("time")
            .astype("float32")
        )
        add_bounds(ds)

        # Add attributes
        for var in PREPROCESSED_DATA_ATTRIBUTES:
            ds[var].attrs.update(PREPROCESSED_DATA_ATTRIBUTES[var])
        ds.attrs.update(input_data_attrs)

        # Save preprocessed data
        filename = f"{datetime.strftime('%Y-%m-%d')}_fluxes_storages.nc"
        output_path = config.preprocessed_data_folder / filename

        if is_new_day:
            comp = dict(zlib=True, complevel=8)
            encoding = {var: comp for var in ds.data_vars}
            time_encoding = {"units": "seconds since 1900-01-01"}
            encoding["time"] = time_encoding
            ds.to_netcdf(
                output_path, unlimited_dims=["time"], mode="w", encoding=encoding
            )

        else:
            append_to_netcdf(output_path, ds, expanding_dim="time")
