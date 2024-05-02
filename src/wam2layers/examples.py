"""Functions for downloading example cases."""

import logging
import shutil
from pathlib import Path
from typing import Callable

import requests
import threddsclient

logger = logging.getLogger(__name__)


def get_4tu_id(doi):
    """Extract 4TU ID from DOI.

    As far as I can tell, DOI is constructed like so:
    10.4121/<4TU-ID>.v<version>

    For lack of a better method, use this to reconstruct 4TU ID
    """
    return doi.split("/")[-1].split(".")[0]


def get_4TU_catalog(uuid):
    """Retrieve the catalogue from 4TU API.

    API is undocumented, but similar to figshare.
    Catalog is available via articles endpoint, under
    custom_fields -> Data Link
    """
    response = requests.get(f"https://data.4tu.nl/v2/articles/{uuid}")
    response.raise_for_status()
    dataset_info = response.json()

    for field in dataset_info["custom_fields"]:
        if field["name"] == "Data Link":
            catalog_url_html = field["value"][0]

    catalog_url_xml = catalog_url_html.replace(".html", ".xml")
    return catalog_url_xml


def create_target_directory(name):
    """Create the target directory but prevent overwrites."""
    target_directory = Path.cwd() / name
    logger.info(f"Creating directory {target_directory}.")

    try:
        target_directory.mkdir(exist_ok=False)
    except FileExistsError as exc:
        raise FileExistsError(
            "Target directory exists. Stopping download to prevent overwriting your files. "
            "If you want to replace your files, first (re)move the directory manually."
        ) from exc

    return target_directory


def download(url, destination):
    """Download file from url to destination."""
    local_filename = url.split("/")[-1]
    local_path = destination / local_filename

    logging.info(f"Downloading {local_filename}")
    with requests.get(url, stream=True) as r:
        with open(local_path, "wb") as f:
            shutil.copyfileobj(r.raw, f)


def download_from_4TU(doi, name):
    """Download each file using threddsclient as crawler."""
    logger.info(f"Downloading 4TU dataset to ./{name}")

    target_dir = create_target_directory(name)
    id_4tu = get_4tu_id(doi)
    catalog_url = get_4TU_catalog(id_4tu)

    # Extract individual file urls from catalogue and download them separately
    file_urls = threddsclient.download_urls(catalog_url)
    for file_url in file_urls:
        download(file_url, target_dir)


AVAILABLE_CASES: dict[str, Callable[[], None]] = {
    "example-volta": "10.4121/bbe10a2a-39dc-4098-a69f-0f677d06ecdd.v1",
    "example-eiffel": "10.4121/f9572240-f179-4338-9e1b-82c5598529e2.v1",
}
