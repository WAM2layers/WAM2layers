"""Prepare test data using ERA5 model level fields.

This script processes daily ERA5 field on multiple model levels.
The data can be downloaded via CDS. More information about the
variables and fields that are needed for tests can be found in
our example ERA5 data downloading script in "scripts/download_era5_ml.py".
"""

from pathlib import Path

import xarray as xr

path_to_data = Path("./path_to_daily_ERA5_data")
era5_name_pattern = "ERA5_2022-08-31"
variable_list = ["e", "sp", "tp", "ml_q", "ml_u", "ml_v"]

# crop raw ERA5 field
for variable in variable_list:
    print("Cropping", f"{era5_name_pattern}_{variable}")
    ds = xr.open_dataset(path_to_data / f"{era5_name_pattern}_{variable}.nc")
    cutout = ds.sel(latitude=slice(55, 45), longitude=slice(0, 10))
    cutout.attrs[
        "Licence"
    ] = "Generated using Copernicus Climate Change Service information [2023]."
    cutout.to_netcdf(f"{era5_name_pattern}_{variable}.nc")

# create source region from the test data
ds_tp = xr.open_dataset(f"{era5_name_pattern}_tp.nc")
source_region = xr.zeros_like(ds_tp["tp"].isel(time=0)).rename("source_region")
source_region[15:35, 15:35] = 1
source_region.to_netcdf("./source_region.nc")
