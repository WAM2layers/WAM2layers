"""Tests for entire workflows."""
from pathlib import Path

import matplotlib.testing.compare
import numpy
import xarray as xr
from click.testing import CliRunner

from wam2layers.cli import cli


def test_preprocess():
    runner = CliRunner()
    result = runner.invoke(
        cli, ["preprocess", "era5", "tests/test_data/test_config.yaml"]
    )
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
    result = runner.invoke(cli, ["backtrack", "tests/test_data/test_config.yaml"])
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


def test_visualize():
    runner = CliRunner()
    result = runner.invoke(
        cli, ["visualize", "output", "tests/test_data/test_config.yaml"]
    )
    assert result.exit_code == 0
    output_path = Path("tests/tmp/output_data/figures/cumulative_sources.png")
    assert output_path.exists()
    expected_output = Path("tests/test_data/verify_output/cumulative_sources.png")
    stdout = matplotlib.testing.compare.compare_images(
        expected_output, output_path, tol=5
    )
    if stdout:
        raise ValueError(
            "Output figures are different! Please check the difference in the produced \n"
            " `cumulative_sources-failed-diff.png` figure and verify your results. \n"
            "If you want to keep the new results and include it in your commit, \n"
            "you can update the reference figure with the following command: \n"
            "`cp tests/tmp/output_data/figures/cumulative_sources.png tests/test_data/verify_output/cumulative_sources.png`"
        )
