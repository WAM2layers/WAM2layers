"""Functions for downloading example cases."""

import logging
import re
import shutil
from pathlib import Path

import requests
import threddsclient

logger = logging.getLogger(__name__)


def parse_doi(doi_string):
    """Extract prefix and suffix from DOI.

    Structure of DOIs:

        <resolver>/<prefix>/<suffix>

    For example:

        https://doi.org/10.4121/f9572240-f179-4338-9e1b-82c5598529e2.v1

    Here the prefix is `10.4121` is specific for 4TU. The suffix (generated by
    4TU) is `f9572240-f179-4338-9e1b-82c5598529e2.v1`

    Note that the suffix may contain slashes, and 4TU inserts a version in it,
    but this is not part of the official doi spec.

    For more info, see
    https://www.crossref.org/documentation/member-setup/constructing-your-dois/
    """
    regex = r"(?P<prefix>10.\d{4,9})\/(?P<suffix>[-._;()\/:A-Z0-9]+$)"
    matches = re.search(regex, doi_string, re.IGNORECASE)

    if matches:
        # Extract the capturing groups
        prefix = matches.group("prefix")
        suffix = matches.group("suffix")

        return prefix, suffix
    else:
        raise ValueError(f"Unable to parse DOI {doi_string}.")


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


def download(url, destination, local_filename=None):
    """Download file from url to destination."""
    if local_filename is None:
        local_filename = url.split("/")[-1]

    local_path = destination / local_filename

    logging.info(f"Downloading {local_filename}")
    with requests.get(url, stream=True) as r:
        with open(local_path, "wb") as f:
            shutil.copyfileobj(r.raw, f)


def get_4tu_id(doi):
    """Extract 4TU ID from DOI.

    As far as I can tell, DOI is constructed like so:
    10.4121/<4TU-ID>.v<version>

    Where 10.4121 is the prefix for 4TU.
    """
    _, suffix = parse_doi(doi)

    # Get UUID from DOI
    id_4tu = suffix.split(".")[0]

    return id_4tu


def get_4tu_catalog(uuid):
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


def get_4tu_files(uuid):
    """Retrieve file names and download links for files.

    This concerns any (additional) files that are not available via OpenDAP
    (e.g. config.yaml).
    """
    response = requests.get(f"https://data.4tu.nl/v2/articles/{uuid}")
    response.raise_for_status()
    dataset_info = response.json()

    files = {}
    for file_info in dataset_info["files"]:
        name = file_info["name"]
        url = file_info["download_url"]
        files[name] = url

    return files


def download_from_4TU(doi, name):
    """Download each file using threddsclient as crawler.

    Args:
        doi: doi of the dataset (with or without https://doi.org/)
        name: human-readable name for the dataset, used to create target directory
    """
    logger.info(f"Downloading 4TU dataset to ./{name}")

    target_dir = create_target_directory(name)
    id_4tu = get_4tu_id(doi)
    catalog_url = get_4tu_catalog(id_4tu)

    # Extract individual file urls from catalogue and download them separately
    file_urls = threddsclient.download_urls(catalog_url)
    for file_url in file_urls:
        download(file_url, target_dir)

    # Extract additional files from file listing
    additional_files = get_4tu_files(id_4tu)
    for filename, url in additional_files.items():
        download(url, target_dir, filename)


def download_from_doi(doi, name):
    """Parse doi and redirect to provider-specific download function.

    Args:
        doi: doi of the dataset (with or without https://doi.org/)
        name: human-readable name for the dataset, used to create target directory
    """
    prefix, _ = parse_doi(doi)

    if prefix == "10.4121":  # 4TU
        download_from_4TU(doi, name)
    elif prefix == "10.5281":  # Zenodo
        "Download utility for Zenodo has not been implemented."
    else:
        raise ValueError(f"Unknown provider with doi prefix {prefix}, cannot download.")
