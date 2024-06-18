import pandas as pd
from wam2layers.config import Config
from wam2layers.preprocessing.era5 import PREPROCESS_ATTRS, get_dp_modellevels, load_data, log_once, logger, preprocess_precip_and_evap
from wam2layers.preprocessing.pressure_levels import extend_pressurelevels, interp_dp_midpoints
from wam2layers.preprocessing.utils import add_bounds, calculate_humidity
from wam2layers.preprocessing.xarray_append import append_to_netcdf
from wam2layers.reference.variables import ERA5_INPUT_DATA_ATTRIBUTES


import xarray as xr


def get_input_dates(config):
    """Dates for pre-processing."""
    return pd.date_range(
        start=config.preprocess_start_date,
        end=config.preprocess_end_date,
        freq=config.input_frequency,
    )


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
    config = Config.from_yaml(config_file)

    for datetime in get_input_dates(config):
        is_new_day = datetime == datetime.floor("1d")
        logger.info(datetime)

        # 4d fields
        levels = config.levels

        input_data = xr.Dataset()
        input_data["q"] = load_data("q", datetime, config)  # in kg kg-1
        input_data["u"] = load_data("u", datetime, config)  # in m/s
        input_data["v"] = load_data("v", datetime, config)  # in m/s

        surface_data = xr.Dataset()
        surface_data["ps"] = load_data("sp", datetime, config)  # in Pa

        if config.level_type == "model_levels":
            dp = get_dp_modellevels(surface_data["ps"], levels)

        if config.level_type == "pressure_levels":
            d2m = load_data("d2m", datetime, config)  # Dew point in K
            surface_data["qs"] = calculate_humidity(d2m, surface_data["ps"])  # kg kg-1
            surface_data["us"] = load_data("u10", datetime, config)  # in m/s
            surface_data["vs"] = load_data("v10", datetime, config)  # in m/s

            input_data, pb = extend_pressurelevels(input_data, surface_data)
            input_data, dp = interp_dp_midpoints(input_data, surface_data["ps"])

        # Calculate column water vapour
        g = 9.80665  # gravitational accelleration [m/s2]
        cwv = input_data["q"] * dp / g  # (kg/m2)

        try:
            # Calculate column water instead of column water vapour
            tcw = load_data("tcw", datetime, config)  # kg/m2
            correction = tcw / cwv.sum(dim="level")
            cw = correction * cwv  # column water (kg/m2)
            if is_new_day:
                logger.info(
                    "Total column water correction:\n"
                    "    ratio total column water / computed integrated water vapour\n"
                    f"    mean over grid for this timestep {correction.mean().item():.4f}"
                )

        except FileNotFoundError:
            # Fluxes will be calculated based on the column water vapour
            cw = cwv
            log_once(
                logger, "Total column water not available; using water vapour only"
            )

        # Integrate fluxes and states to upper and lower layer
        if config.level_type == "model_levels":
            # TODO: Check if this is a reasonable choice for boundary
            boundary = 111
            lower_layer = dp["level"] > boundary
            upper_layer = ~lower_layer

        if config.level_type == "pressure_levels":
            upper_layer = input_data["p"] < pb.broadcast_like(input_data["p"])
            lower_layer = ~upper_layer

        # Vertically integrate state over two layers
        s_lower = cw.where(lower_layer).sum(dim="level")
        s_upper = cw.where(upper_layer).sum(dim="level")

        # Determine the fluxes
        fx = input_data["u"] * cw  # eastward atmospheric moisture flux (kg m-1 s-1)
        fy = input_data["v"] * cw  # northward atmospheric moisture flux (kg m-1 s-1)

        # Vertically integrate fluxes over two layers
        fx_lower = fx.where(lower_layer).sum(dim="level")  # kg m-1 s-1
        fy_lower = fy.where(lower_layer).sum(dim="level")  # kg m-1 s-1
        fx_upper = fx.where(upper_layer).sum(dim="level")  # kg m-1 s-1
        fy_upper = fy.where(upper_layer).sum(dim="level")  # kg m-1 s-1

        precip, evap = preprocess_precip_and_evap(datetime, config)

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
                    "evap": evap,
                    "precip": precip,
                }
            )
            .expand_dims("time")
            .astype("float32")
        )

        # Verify that the data meets all the requirements for the model
        # check_input(ds)

        add_bounds(ds)

        # Add attributes
        for var in ERA5_INPUT_DATA_ATTRIBUTES:
            ds[var].attrs.update(ERA5_INPUT_DATA_ATTRIBUTES[var])
        ds.attrs.update(PREPROCESS_ATTRS)

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