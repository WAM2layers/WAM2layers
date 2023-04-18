from pathlib import Path

from click.testing import CliRunner

from wam2layers.cli import cli


def test_preprocess():
    runner = CliRunner()
    result = runner.invoke(cli, 'preprocess era5 sample_data/sample_config.yaml')
    assert result.exit_code == 0
    assert len(result.output)
    assert Path("/tmp/wam2layers/sample/preprocessed_data/2021-07-15_fluxes_storages.nc").exists()


def test_backtrack():
    runner = CliRunner()
    result = runner.invoke(cli, 'backtrack sample_data/sample_config.yaml')
    assert result.exit_code == 0
    assert "Experiment complete" in result.output


# TODO: add tests to verify that the data hasn't changed
