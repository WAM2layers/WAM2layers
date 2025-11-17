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
from pathlib import Path

from wam2layers.config import Config
from wam2layers.tracking.io import input_path

# use HTTP endpoint with curl? But can't subset in that case.
OPENDAP_ENDPOINT_HTTP = "https://opendap.4tu.nl/thredds/fileServer/data2/djht/00f7fa45-899e-4573-ae23-234f6c5193d0/1"

# TODO: make it work with DAP endpoint so we can do the spatial subset. However
# requires complex subsetting that is exactly why it would be so much easier
# with xarray. And then we might as well do it on the fly.
OPENDAP_ENDPOINT_DAP4 = "dap4://opendap.4tu.nl/thredds/dap4/data2/djht/00f7fa45-899e-4573-ae23-234f6c5193d0/1"


def download_file_http(url: str, output_path: Path):
    """Download a file using curl."""
    cmd = ["curl", "-L", "-o", str(output_path), url]
    print(f"Downloading {url} -> {output_path}")
    subprocess.run(cmd, check=True)


def read_file_with_opendap(url):
    """Download a file using opendap workaround.

    Uses pydap with DAP4 backend seems to work okayish with recent version of xarray
    perhaps due to https://github.com/pydata/xarray/pull/10482 ??
    """
    import xarray as xr

    from wam2layers.tracking.io import select_subdomain

    ds = xr.open_dataset(url, engine="netcdf4", decode_cf=False)

    # 2. Drop string vars to avoid segfault # TODO: check if necessary
    string_vars = [v for v in ds.data_vars if ds[v].dtype.kind in ("O", "S", "U")]
    ds = ds.drop_vars(string_vars)

    # 3. Sanitize attributes
    for attr in list(ds.attrs):
        if attr.startswith("_dap4") or attr == "_NCProperties":
            del ds.attrs[attr]

    # TODO insert bounding box
    bbox = ...
    ds = select_subdomain(ds, bbox)

    # TODO update filename
    ds.to_netcdf("local_copy.nc")


def download_from_config(config_file: str):
    # Load the WAM2layers config
    cfg = Config.from_yaml(config_file)

    # Extract start/end dates
    start_date = cfg.preprocess_start_date
    end_date = cfg.preprocess_end_date

    # Get the tracking region if needed
    tracking_domain = cfg.tracking_domain
    print(f"Tracking region: {tracking_domain}")

    # Directory to save files
    output_dir = cfg.preprocessed_data_folder
    output_dir.mkdir(parents=True, exist_ok=True)

    # Generate list of dates
    current_date = start_date
    while current_date <= end_date:
        # Construct the URL
        url = input_path(current_date, OPENDAP_ENDPOINT_HTTP)
        output_file = output_dir / url.split("/")[-1]  # save with filename only
        download_file_http(url, output_file)
        current_date += timedelta(days=1)


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        print("Usage: python download_wam2layers.py <config.yaml>")
        sys.exit(1)

    main(sys.argv[1])
