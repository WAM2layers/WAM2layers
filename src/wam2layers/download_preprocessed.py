"""Download preprocessed WAM2layers input data

This module provides functionality to download pre-processed data from 4TU based
on a WAM2layers config file. It will read the input date range, and extent of
the tracking domain, and download only the necessary data.

It can be invoked through

    wam2layers download-preprocessed example-config.yaml

where you replace "example-config.yaml" with your own configuration.
"""

import subprocess
from datetime import timedelta

import xarray as xr

from wam2layers.config import Config
from wam2layers.tracking.io import select_subdomain

ENDPOINT_HTTP = "https://opendap.4tu.nl/thredds/fileServer/data2/djht/00f7fa45-899e-4573-ae23-234f6c5193d0/1"
ENDPOINT_DAP4 = "dap4://opendap.4tu.nl/thredds/dap4/data2/djht/00f7fa45-899e-4573-ae23-234f6c5193d0/1"


def get_filename(date):
    return f"{date.strftime('%Y-%m-%d')}_fluxes_storages.nc"


def download_http(date, output_dir, bbox: str | None = None):
    """Download a file using curl.

    This downloads the entire domain, then crops the domain.
    """
    filename = get_filename(date)
    url = f"{ENDPOINT_HTTP}/{filename}"
    output_file = output_dir / filename

    cmd = ["curl", "-L", "-o", str(output_file), url]
    subprocess.run(cmd, check=True)

    if bbox is not None:
        print("Cropping tracking domain...", end=" ", flush=True)
        with xr.open_dataset(output_file) as ds_full:
            # This automatically closes the file after reading the data
            ds_cropped = select_subdomain(ds_full, bbox).load()

        # Overwrite the original output file
        ds_cropped.to_netcdf(output_file)
        print("Done")


def download_dap4(date, output_dir, bbox: str | None = None):
    """Download a file using opendap workaround.

    This can subset the domain on the server, resulting in smaller data
    transfers, but the protocol is less stable.
    """
    # DAP4 works now, but returns garbage data.
    # Specifically, e.g. looking at ds.precip.values, initial chunks look good,
    # but later ones filled with rubbish.
    raise UserWarning(
        "DAP4 download seems unreliable and should not be used at the moment."
    )

    filename = get_filename(date)
    url = f"{ENDPOINT_DAP4}/{filename}"
    output_file = output_dir / filename

    ds = xr.open_dataset(url, engine="pydap", decode_times=False)

    # Sanitize attributes
    for attr in list(ds.attrs):
        if attr.startswith("_dap4") or attr == "_NCProperties":
            del ds.attrs[attr]

    if bbox:
        ds = select_subdomain(ds, bbox)

    ds.to_netcdf(output_file)


def download_from_config(config_file: str, protocol="http", crop=True):
    """Download pre-processed data from 4TU server to local directory.

    - reads dates from input config file
    - downloads the corresponding data
    - crops the data to the tracking domain (optional)
    - stores data in `preprocessed_data_path`

    Choose between `http` and `dap4` backends.
    - http is more stable but downloads the full file and (optionally) crops afterwards.
    - dap4 can subset on the server, saving on data transfers, but is less stable.
    """
    # Load the WAM2layers config
    cfg = Config.from_yaml(config_file)

    # Extract start/end dates
    start_date = cfg.preprocess_start_date
    end_date = cfg.preprocess_end_date

    # Directory to save files
    output_dir = cfg.preprocessed_data_folder
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Downloading WAM2layers data to {output_dir}")

    # Get the tracking domain if needed
    if crop:
        bbox = str(cfg.tracking_domain)
        print(f"Tracking domain: {bbox}")
        if bbox == "None":
            bbox = None
    else:
        bbox = None

    current_date = start_date
    while current_date <= end_date:
        print(f"Downloading data for {current_date}")

        if protocol == "http":
            download_http(current_date, output_dir, bbox)
        else:  # protocol == "dap4"
            download_dap4(current_date, output_dir, bbox)

        current_date += timedelta(days=1)
