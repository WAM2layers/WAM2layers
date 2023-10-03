"""Tests for entire workflows."""
from pathlib import Path
from click.testing import CliRunner

from wam2layers.cli import cli


def test_preprocess():
    runner = CliRunner()
    result = runner.invoke(cli, ["preprocess", "era5", "tests/test_config.yaml"])
    assert result.exit_code == 0
    assert Path("tests/tmp/preprocessed_data/2022-08-31_fluxes_storages.nc").exists()
    # TODO: verify output


def test_backtrack():
    runner = CliRunner()
    result = runner.invoke(cli, ["backtrack", "tests/test_config.yaml"])
    assert result.exit_code == 0
    assert Path("tests/tmp/output_data/backtrack_2022-08-31T18-00.nc").exists()
    # TODO: verify output
