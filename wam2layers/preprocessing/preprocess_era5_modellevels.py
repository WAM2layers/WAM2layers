import os

import pandas as pd
import xarray as xr
import yaml
import numpy as np
from pathlib import Path

from preprocessing.preprocessing import get_grid_info

# Read case configuration
with open("../cases/era5_2021.yaml") as f:
    config = yaml.safe_load(f)

# Create the preprocessed data folder if it does not exist yet
output_dir = Path(config["preprocessed_data_folder"]).expanduser()
output_dir.mkdir(exist_ok=True, parents=True)

# Set constants
g = 9.80665  # [m/s2]
density_water = 1000  # [kg/m3]

# Preselection of model levels
modellevels = config["modellevels"]
print("Number of model levels:", len(modellevels))

# Load a and b coefficients
df = pd.read_csv("tableERA5model_to_pressure.csv")
a = df["a [Pa]"].to_xarray().rename(index="lev")  # .sel(lev=modellevels)
b = df["b"].to_xarray().rename(index="lev")  # .sel(lev=modellevels)

# Calculate a and b at mid levels (model levels)
a_full = ((a[1:] + a[:-1].values) / 2.0).sel(lev=modellevels)
b_full = ((b[1:] + b[:-1].values) / 2.0).sel(lev=modellevels)

# Construct a and b at edges for selected levels
a_edge = xr.concat([a[0], (a_full[1:].values + a_full[:-1]) / 2.0, a[-1]], dim="lev")
b_edge = xr.concat([b[0], (b_full[1:].values + b_full[:-1]) / 2.0, b[-1]], dim="lev")


def load_surface_data(variable, date):
    """Load data for given variable and date."""
    filename = f"FloodCase_202107_{variable}.nc"
    filepath = os.path.join(config["input_folder"], filename)
    da = xr.open_dataset(filepath)[variable]

    # Include midnight of the next day (if available)
    extra = date + pd.Timedelta(days=1)
    return da.sel(time=slice(date, extra))


def load_modellevel_data(variable, date):
    """Load model level data for given variable and date."""
    filename = f"FloodCase_202107_ml_{variable}.nc"
    filepath = os.path.join(config["input_folder"], filename)
    da = xr.open_dataset(filepath)[variable]

    # Include midnight of the next day (if available)
    extra = date + pd.Timedelta(days=1)
    return da.sel(time=slice(date, extra)).sel(lev=modellevels)


datelist = pd.date_range(
    start=config["preprocess_start_date"],
    end=config["preprocess_end_date"],
    freq="d",
    inclusive="left",
)

for date in datelist[:]:
    print(date)

    # Load data
    u = load_modellevel_data("u", date)
    v = load_modellevel_data("v", date)
    q = load_modellevel_data("q", date)
    sp = load_surface_data("sp", date)  # in Pa
    evap = load_surface_data("e", date)
    cp = load_surface_data("cp", date)
    lsp = load_surface_data("lsp", date)
    precip = cp + lsp

    # Get grid info
    a_gridcell, l_ew_gridcell, l_mid_gridcell = get_grid_info(u)

    # Transfer negative (originally positive) values of evap to precip
    precip = np.maximum(precip, 0) + np.maximum(evap, 0)
    # Change sign convention to all positive,
    evap = np.abs(np.minimum(evap, 0))

    # Convert model levels to pressure values
    sp_modellevels2 = sp.expand_dims({"lev": modellevels}, axis=1)
    p_modellevels = a_edge + b_edge * sp_modellevels2  # in Pa

    # calculate the difference between the pressure levels
    dp_modellevels = p_modellevels.diff(dim="lev")  # in Pa

    # Split in 2 layers
    # To do: Check if this is a reasonable choice
    boundary = 111  # set boundary model level 111

    # Integrate fluxes and states to upper and lower layer
    lower_layer = dp_modellevels.lev > boundary
    upper_layer = ~lower_layer

    cwv = q * dp_modellevels / g  #  # column water vapor (kg/m2)

    if config["vertical_integral_available"] == True:
        # calculate column water instead of column water vapour
        tcw = load_surface_data("tcw", date)  # kg/m2
        cw = (tcw / cwv.sum(dim="lev")) * cwv  # column water (kg/m2)
    else:
        # calculate the fluxes based on the column water vapour
        cw = cwv

    # Vertically integrate state over two layers
    s_lower = (
        cw.where(lower_layer).sum(dim="lev")
        * a_gridcell[None, :, None]
        / density_water
    )  # m3
    s_upper = (
        cw.where(upper_layer).sum(dim="lev")
        * a_gridcell[None, :, None]
        / density_water
    )  # m3

    # Determine the fluxes
    fx = u * cw  # eastward atmospheric moisture flux (kg m-1 s-1)
    fy = v * cw  # northward atmospheric moisture flux (kg m-1 s-1)

    # Vertically integrate fluxes over two layers
    fx_lower = fx.where(lower_layer).sum(dim="lev")  # kg m-1 s-1
    fy_lower = fy.where(lower_layer).sum(dim="lev")  # kg m-1 s-1

    fx_upper = fx.where(upper_layer).sum(dim="lev")  # kg m-1 s-1
    fy_upper = fy.where(upper_layer).sum(dim="lev")  # kg m-1 s-1

    # Save preprocessed data
    filename = f"{date.strftime('%Y-%m-%d')}_fluxes_storages.nc"
    output_path = os.path.join(config["preprocessed_data_folder"], filename)
    xr.Dataset(
        {  # TODO: would be nice to add coordinates and units as well
            "fx_upper": fx_upper,
            "fy_upper": fy_upper,
            "fx_lower": fx_lower,
            "fy_lower": fy_lower,
            "s_upper": s_upper,
            "s_lower": s_lower,
            "evap": evap,
            "precip": precip,
        }
    ).to_netcdf(output_path)
