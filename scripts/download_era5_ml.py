"""Download ERA5.py

This script can be used to download ERA5 data on model levels using the CDS API.
Please see the installation instructions: https://cds.climate.copernicus.eu/api-how-to
You need to have a valid CDS API key and you need to pip install cdsapi.

For model levels the following variables are downloaded:
  - u, v, q on selected levels
  - tp, e, sp, and tcw at the surface

Modify the settings below, then run with:

python download_era5_ml.py
"""
import cdsapi
import pandas as pd
from pathlib import Path

target_dir = "."
skip_exist = True

datelist = pd.date_range("20200101", "20200103")

area = None  # None for global, or [N, W, S, E]
grid = [0.25, 0.25]

times = [
    '00:00',
    '01:00',
    '02:00',
    '03:00',
    '04:00',
    '05:00',
    '06:00',
    '07:00',
    '08:00',
    '09:00',
    '10:00',
    '11:00',
    '12:00',
    '13:00',
    '14:00',
    '15:00',
    '16:00',
    '17:00',
    '18:00',
    '19:00',
    '20:00',
    '21:00',
    '22:00',
    '23:00',
]

levels = list(range(137))

ml_variables = {
    "u": 131,
    "v": 132,
    "q": 133,
}

surface_variables = {
    'tp': "total_precipitation",
    'e': "evaporation",
    'sp': "mean_sea_level_pressure",
    'tcw': "total_column_water",
}

## The part below should not have to be modified
################################################

c = cdsapi.Client()

# We want one file per variable per day
for date in datelist:

    # Create data directory if it doesn't exist yet
    outfolder = Path(target_dir) / str(date.year) / str(date.month)
    outfolder.mkdir(exist_ok=True, parents=True)

    # Download surface variables
    for variable, long_name in surface_variables.items():
        outfile = f"ERA5_{date.strftime('%Y-%m-%d')}_{variable}.nc"
        if (outfolder / outfile).exists() and skip_exist:
            print(
                f"{outfolder / outfile} already exists, skipping. Set skip_exist = False to force re-download"
            )
        else:
            try:
                c.retrieve(
                    'reanalysis-era5-single-levels', {
                        'product_type': 'reanalysis',
                        'variable': long_name,
                        'date': date.strftime("%Y-%m-%d"),
                        'time': times,
                        'area': area,
                        'grid': grid,
                        'format': 'netcdf',
                    }, str(outfolder / outfile))
            except Exception as e:
                print(e)
                print(
                    ""
                    f"Request failed for {variable}, {date}. Proceeding. You can"
                    " run the code again and set 'skip_exist' to true to avoid"
                    " duplicate downloads."
                    "")

    # Download 3d variables
    for variable, param in ml_variables.items():
        outfile = f"ERA5_{date.strftime('%Y-%m-%d')}_{variable}_ml.nc"
        if (outfolder / outfile).exists() and skip_exist:
            print(
                f"{outfolder / outfile} already exists, skipping. Set skip_exist = False to force re-download"
            )
        else:
            try:
                c.retrieve(
                    'reanalysis-era5-complete', {
                        'time': times,
                        'date': date.strftime("%Y-%m-%d"),
                        'levelist': levels,
                        'param': param,
                        'class': 'ea',
                        'expver': '1',
                        'levtype': 'ml',
                        'stream': 'oper',
                        'type': 'an',
                        'format': 'netcdf',
                        'area': area,
                        'grid': grid,
                    }, str(outfolder / outfile))
            except Exception as e:
                print(e)
                print(
                    ""
                    f"Request failed for {variable}, {date}. Proceeding. You can"
                    " run the code again and set 'skip_exist' to true to avoid"
                    " duplicate downloads."
                    "")
