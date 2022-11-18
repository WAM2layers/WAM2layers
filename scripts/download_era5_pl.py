"""Download ERA5.py

This script can be used to download ERA5 data on pressure levels using the CDS API.
Please see the installation instructions: https://cds.climate.copernicus.eu/api-how-to
You need to have a valid CDS API key and you need to pip install cdsapi.

For pressure levels the following variables are downloaded:
  - u, v, q on selected levels
  - tp, e, sp, t2m, d2m, u10, v10, tcw at the surface

Modify the settings below, then run with:

python download_era5_pl.py
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

levels = [
    '100',
    '200',
    '300',
    '400',
    '500',
    '550',
    '600',
    '650',
    '700',
    '750',
    '775',
    '800',
    '825',
    '850',
    '875',
    '900',
    '925',
    '950',
    '975',
    '1000',
]

pl_variables = {
    'q': 'specific_humidity',
    'u': 'u_component_of_wind',
    'v': 'v_component_of_wind',
}

surface_variables = {
    'tp': 'total_precipitation',
    'e': 'evaporation',
    'sp': 'surface_pressure',
    'd2m': '2m_dewpoint_temperature',
    'u10': '10m_u_component_of_wind',
    'v10': '10m_v_component_of_wind',
    'tcw': 'total_column_water',
}

## The part below should not have to be modified
################################################

c = cdsapi.Client()

# We want one file per variable per day
for date in datelist:

    # Create data directory if it doesn't exist yet
    outfolder = Path(target_dir) / f"{date.year:04d}" / f"{date.month:02d}" / f"{date.day:02d}"
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
    for variable, long_name in pl_variables.items():
        outfile = f"ERA5_{date.strftime('%Y-%m-%d')}_pl_{variable}.nc"
        if (outfolder / outfile).exists() and skip_exist:
            print(
                f"{outfolder / outfile} already exists, skipping. Set skip_exist = False to force re-download"
            )
        else:
            try:
                c.retrieve(
                    'reanalysis-era5-pressure-levels', {
                        'time': times,
                        'date': date.strftime("%Y-%m-%d"),
                        'pressure_level': levels,
                        'variable': long_name,
                        'area': area,
                        'grid': grid,
                        'product_type': 'reanalysis',
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
