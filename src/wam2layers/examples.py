"""Functions for preparing example cases.


Each prepare function should:

* Create a new directory in the current working directory
* Add a configuration file (may be shipped with the repo or fetched from elsewhere)
* Add meteorological input data (download from external repository)
* Add a tracking_region file (may be shipped with the repo or fetched from elsewhere)

Sources should be persistent (ideally with DOI)
"""

import io
import zipfile
from pathlib import Path
from typing import Callable

import requests


def extract_from_4TU_archive(download_link, target_dir):
    # See https://stackoverflow.com/a/14260592
    response = requests.get(download_link, stream=True)  # TODO stream necessary or not?
    assert response.ok, "Download failed"
    z = zipfile.ZipFile(
        io.BytesIO(response.content)
    )  # TODO should I use context managers here?
    z.extractall(target_dir)


def download_volta() -> None:
    # TODO: replace with correct link; this is a random other dataset that's already published
    download_link = "https://data.4tu.nl/ndownloader/items/4e291b8f-a37e-4378-8ca6-954a44fdc8fb/versions/2"
    target_directory = Path.cwd() / "example_volta"

    print(
        f"Downloading and extracing example data for Volta forward tracking case to {target_directory}"
    )
    extract_from_4TU_archive(download_link, target_directory)


def download_eiffel() -> None:
    # TODO: replace with correct link; this is a random other dataset that's already published
    download_link = "https://data.4tu.nl/ndownloader/items/4e291b8f-a37e-4378-8ca6-954a44fdc8fb/versions/2"
    target_directory = Path.cwd() / "example_eiffel"

    print(
        f"Downloading and extracing example data for Eiffel backward tracking case to {target_directory}"
    )

    extract_from_4TU_archive(download_link, target_directory)


AVAILABLE_CASES: dict[str, Callable[[], None]] = {
    "example-volta": download_volta,
    "example-eiffel": download_eiffel,
}
