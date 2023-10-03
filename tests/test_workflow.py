"""Tests for entire workflows."""
import numpy
from pathlib import Path
import xarray as xr
from click.testing import CliRunner

from wam2layers.cli import cli


def test_preprocess():
    runner = CliRunner()
    result = runner.invoke(cli, ["preprocess", "era5", "tests/test_config.yaml"])
    assert result.exit_code == 0
    output_path = Path("tests/tmp/preprocessed_data/2022-08-31_fluxes_storages.nc")
    assert output_path.exists()
    # verify outputs
    expected_output = xr.open_dataset(
        "tests/test_data/verify_output/2022-08-31_fluxes_storages.nc"
    )
    output = xr.open_dataset(output_path)
    numpy.testing.assert_allclose(
        expected_output["precip"].values, output["precip"].values
    )


def test_backtrack():
    runner = CliRunner()
    result = runner.invoke(cli, ["backtrack", "tests/test_config.yaml"])
    assert result.exit_code == 0
    output_path = Path("tests/tmp/output_data/backtrack_2022-08-31T18-00.nc")
    assert output_path.exists()
    # verify outputs
    expected_output = xr.open_dataset(
        "tests/test_data/verify_output/backtrack_2022-08-31T18-00.nc"
    )
    output = xr.open_dataset(output_path)
    numpy.testing.assert_allclose(
        expected_output["e_track"].values, output["e_track"].values
    )
