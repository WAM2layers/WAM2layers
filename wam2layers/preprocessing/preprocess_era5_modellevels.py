from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import yaml

from checks import check_input  # TODO move to other module

# Set constants
g = 9.80665  # [m/s2]

# Read case configuration
with open("../../cases/era5_2021.yaml") as f:
    config = yaml.safe_load(f)

# Create the preprocessed data folder if it does not exist yet
output_dir = Path(config["preprocessed_data_folder"]).expanduser()
output_dir.mkdir(exist_ok=True, parents=True)


def load_data(variable, date, levels=None):
    """Load data for given variable and date."""
    # TODO: remove hardcoded filename pattern
    prefix = "FloodCase_202107_ml" if levels is not None else "FloodCase_202107"
    filepath = Path(config["input_folder"]) / f"{prefix}_{variable}.nc"
    da = xr.open_dataset(filepath)[variable]

    # Include midnight of the next day (if available)
    extra = date + pd.Timedelta(days=1)

    if levels is not None:
        return da.sel(time=slice(date, extra)).sel(lev=levels)
    return da.sel(time=slice(date, extra))


# Preselection of model levels
modellevels = config["modellevels"]
print("Number of model levels:", len(modellevels))

# Load a and b coefficients
df = pd.read_csv("tableERA5model_to_pressure.csv")
a = df["a [Pa]"].to_xarray().rename(index="lev")
b = df["b"].to_xarray().rename(index="lev")

# Calculate a and b at mid levels (model levels)
a_full = ((a[1:] + a[:-1].values) / 2.0).sel(lev=modellevels)
b_full = ((b[1:] + b[:-1].values) / 2.0).sel(lev=modellevels)

# Construct a and b at edges for selected levels
a_edge = xr.concat([a[0], (a_full[1:].values + a_full[:-1]) / 2, a[-1]], dim="lev")
b_edge = xr.concat([b[0], (b_full[1:].values + b_full[:-1]) / 2, b[-1]], dim="lev")

datelist = pd.date_range(
    start=config["preprocess_start_date"],
    end=config["preprocess_end_date"],
    freq="d",
    inclusive="left",
)

for date in datelist[:]:
    print(date)

    # 4d fields
    u = load_data("u", date, levels=modellevels)  # in m/s
    v = load_data("v", date, levels=modellevels)  # in m/s
    q = load_data("q", date, levels=modellevels)  # in kg kg-1

    # Precipitation and evaporation
    evap = load_data("e", date)  # in m (accumulated hourly)
    cp = load_data("cp", date)  # convective precip in m (accumulated hourly)
    lsp = load_data("lsp", date)  # large scale precip in m (accumulated)
    precip = (cp + lsp)

    # 3d fields
    p_surf = load_data("sp", date)  # in Pa

    # Transfer negative (originally positive) values of evap to precip
    precip = np.maximum(precip, 0) + np.maximum(evap, 0)
    precip = precip.assign_attrs(units="m")  # xr doesn't preserve attributes in arithmetic
    # Change sign convention to all positive,
    evap = np.abs(np.minimum(evap, 0))

    # Split in 2 layers
    # Convert model levels to pressure values
    sp_modellevels2 = p_surf.expand_dims({"lev": modellevels}, axis=1)
    p_modellevels = a_edge + b_edge * sp_modellevels2  # in Pa

    # calculate the difference between the pressure levels
    dp_modellevels = p_modellevels.diff(dim="lev")  # in Pa
    # To do: Check if this is a reasonable choice
    boundary = 111  # set boundary model level 111

    cwv = q * dp_modellevels / g  # column water vapor (kg/m2)

    if config["vertical_integral_available"]:
        # calculate column water instead of column water vapour
        tcw = load_data("tcw", date)  # kg/m2
        cw = (tcw / cwv.sum(dim="lev")) * cwv  # column water (kg/m2)
        # TODO: warning if cw >> cwv
    else:
        # fluxes will be calculated based on the column water vapour
        cw = cwv
    cw = cw.assign_attrs(units="kg m-2")

    # Integrate fluxes and states to upper and lower layer
    lower_layer = dp_modellevels.lev > boundary
    upper_layer = ~lower_layer

    # Vertically integrate state over two layers
    s_lower = cw.where(lower_layer).sum(dim="lev", keep_attrs=True)
    s_upper = cw.where(upper_layer).sum(dim="lev", keep_attrs=True)

    # Determine the fluxes
    fx = (u * cw).assign_attrs(units="kg m-1 s-1")  # eastward atmospheric moisture flux (kg m-1 s-1)
    fy = (v * cw).assign_attrs(units="kg m-1 s-1")  # northward atmospheric moisture flux (kg m-1 s-1)

    # Vertically integrate fluxes over two layers
    fx_lower = fx.where(lower_layer).sum(dim="lev", keep_attrs=True)  # kg m-1 s-1
    fy_lower = fy.where(lower_layer).sum(dim="lev", keep_attrs=True)  # kg m-1 s-1
    fx_upper = fx.where(upper_layer).sum(dim="lev", keep_attrs=True)  # kg m-1 s-1
    fy_upper = fy.where(upper_layer).sum(dim="lev", keep_attrs=True)  # kg m-1 s-1

    # Combine everything into one dataset
    ds = xr.Dataset(
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
    )

    # Verify that the data meets all the requirements for the model
    check_input(ds)

    # Save preprocessed data
    filename = f"{date.strftime('%Y-%m-%d')}_fluxes_storages.nc"
    output_path = output_dir / filename
    ds.to_netcdf(output_path)
