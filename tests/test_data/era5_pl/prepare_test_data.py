"""Prepare test data using ERA5 model level fields.

This script processes daily ERA5 field on multiple model levels.
The data can be downloaded via CDS. More information about the
variables and fields that are needed for tests can be found in
our example ERA5 data downloading script in "scripts/download_era5_ml.py".
"""

from pathlib import Path

import xarray as xr

path_to_data = Path("raw/")
era5_name_pattern = "ERA5_2020-01-01"
variable_list = ["d2m", "e", "pl_q", "pl_u", "pl_v", "sp", "tcw", "tp", "u10", "v10"]

# crop raw ERA5 field
for variable in variable_list:
    print("Cropping", f"{era5_name_pattern}_{variable}")
    ds = xr.open_dataset(path_to_data / f"{era5_name_pattern}_{variable}.nc")
    cutout = ds.sel(latitude=slice(55, 45), longitude=slice(0, 10))
    cutout.attrs[
        "Licence"
    ] = "Generated using Copernicus Climate Change Service information [2023]."
    encoding = {var: ds[var].encoding.update(zlib=True) for var in ds.data_vars}
    cutout.to_netcdf(f"{era5_name_pattern}_{variable}.nc")

