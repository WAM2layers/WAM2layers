"""Preprocessing input data validation."""
import numpy as np
import xarray as xr

from wam2layers.config import Config


class ValidationError(Exception):
    ...


def validate_input(data: xr.Dataset, config: Config) -> None:
    validate_dims(data)
    validate_vars(data, config)
    check_nans(
        data
    )  # note: needs to load all data in memory... move to pre-output writing


def validate_dims(data: xr.Dataset) -> None:
    expected_dims = ["latitude", "longitude", "level"]
    for dim in expected_dims:
        if dim not in data.dims:
            msg = (
                f"Dimension '{dim}' is missing from the input data.\n"
                f" Dims present in data: {data.dims}."
            )
            raise ValidationError(msg)


def validate_vars(data: xr.Dataset, config: Config) -> None:
    expected_level_vars = ["q", "u", "v"]
    expected_surface_vars = ["ps", "precip", "evap"]
    if config.level_type == "pressure_levels":
        expected_surface_vars.extend(["qs", "us", "vs"])
    if config.level_type == "model_levels":
        expected_level_vars.extend(["dp"])

    # Check if the variables are present
    for var in expected_level_vars + expected_surface_vars:
        if var not in data.data_vars:
            msg = (
                f"Variable '{var}' is missing from the input data.\n"
                f" Variables present in data: {data.data_vars}."
            )
            raise ValidationError(msg)

    # Make sure surface variables have correct dims
    for var in expected_surface_vars:
        if set(data[var].dims) != {"latitude", "longitude"}:
            msg = (
                "Surface variables should only have the dimensions latitude, "
                "longitude.\n"
                f" Var {var} has dimensions: {data[var].dims}."
            )
            raise ValidationError(msg)

    for var in expected_level_vars:
        if set(data[var].dims) != {"level", "latitude", "longitude"}:
            msg = (
                "Pressure/model Level variables should have the dimensions level, "
                "latitude, longitude.\n"
                f" Var {var} has dimensions: {data[var].dims}."
            )
            raise ValidationError(msg)


def check_nans(data: xr.Dataset):
    for var in data.data_vars:
        nans = data[var].isnull()
        if nans.any():
            msg = f"{np.sum(nans)} NaN values found in variable '{var}'"
            raise ValidationError(msg)
