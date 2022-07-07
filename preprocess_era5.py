import os
import numpy as np
import pandas as pd
import xarray as xr
import yaml

from preprocessing import get_grid_info


# Set constants
g = 9.80665  # [m/s2]
density_water = 1000  # [kg/m3]

# Read case configuration
with open("cases/era5_2013_local.yaml") as f:
    config = yaml.safe_load(f)


def load_data(variable, date):
    """Load data for given variable and date."""
    filename = f"FloodCase_201305_{variable}.nc"
    filepath = os.path.join(config["input_folder"], filename)
    da = xr.open_dataset(filepath)[variable]

    # Include midnight of the next day (if available)
    extra = date + pd.Timedelta(days=1)
    return da.sel(time=slice(date, extra))


datelist = pd.date_range(
    start=config["preprocess_start_date"],
    end=config["preprocess_end_date"],
    freq="d",
    inclusive="left",
)

for date in datelist[:]:
    print(date)

    # Load data
    u = load_data("u", date)  # m/s
    v = load_data("v", date)  # m/s
    q = load_data("q", date)  # kg kg-1
    sp = load_data("sp", date)  # Pa
    tcw = load_data("tcw", date) # kg/m2

    # Get grid cell area
    a_gridcell, _, _ = get_grid_info(u)

    # Load evap & precip
    evap = load_data("e", date)  # m (accumulated hourly)
    cp = load_data("cp", date)  # convective precipitation in m (accumulated hourly)
    lsp = load_data("lsp", date)  # large scale precipitation in m (accumulated hourly)
    precip = cp + lsp

    # Transfer negative (originally positive) values of evap to precip
    precip = np.maximum(precip, 0) + np.maximum(evap, 0)
    # Change sign convention to all positive,
    evap = np.abs(np.minimum(evap, 0))

    # Create pressure array
    levels = q.level #in hPa
    p = levels.broadcast_like(u) * 100  # Pa

    # Interpolate to new levels
    edges = 0.5 * (levels.values[1:] + levels.values[:-1])
    u = u.interp(level=edges)
    v = v.interp(level=edges)
    q = q.interp(level=edges)

    # Calculate pressure jump
    dp = p.diff(dim="level")
    dp["level"] = edges

    # Determine the fluxes and states
    fx = u * q * dp / g  # eastward atmospheric moisture flux (kg*m-1*s-1)
    fy = v * q * dp / g  # northward atmospheric moisture flux (#kg*m-1*s-1)
    cwv = q * dp / g * a_gridcell / density_water  # column water vapor (m3)

    # Split in 2 layers
    P_boundary = 0.72878581 * sp + 7438.803223
    lower_layer = (dp.level < sp / 100) & (dp.level > P_boundary / 100)
    upper_layer = dp.level < P_boundary / 100

    # Integrate fluxes and state
    fx_lower = fx.where(lower_layer).sum(dim="level")  # kg*m-1*s-1
    fy_lower = fy.where(lower_layer).sum(dim="level")  # kg*m-1*s-1
    s_lower = cwv.where(lower_layer).sum(dim="level")  # m3

    fx_upper = fx.where(upper_layer).sum(dim="level")  # kg*m-1*s-1
    fy_upper = fy.where(upper_layer).sum(dim="level")  # kg*m-1*s-1
    s_upper = cwv.where(upper_layer).sum(dim="level")  # m3

    print(
        "Check calculation water vapor, this value should be zero:",
        (cwv.sum(dim="level") - (s_upper + s_lower)).sum().values,
    )

    tcwm3 = tcw * a_gridcell[np.newaxis, :] / density_water  # m3

    # Save preprocessed data as daily files including the following midnight
    filename = f"{date.strftime('%Y-%m-%d')}_fluxes_storages.nc"
    output_path = os.path.join(config["preprocessed_data_folder"], filename)
    xr.Dataset({
        "fx_lower": fx_lower.assign_attrs(units="kg*m-1*s-1"),
        "fy_lower": fy_lower.assign_attrs(units="kg*m-1*s-1"),
        "fx_upper": fx_upper.assign_attrs(units="kg*m-1*s-1"),
        "fy_upper": fy_upper.assign_attrs(units="kg*m-1*s-1"),
        "s_upper": s_upper.assign_attrs(units="m**3"),
        "s_lower": s_lower.assign_attrs(units="m**3"),
        "precip": precip.assign_attrs(units="m"),
        "evap": evap.assign_attrs(units="m"),
    }).to_netcdf(output_path)
