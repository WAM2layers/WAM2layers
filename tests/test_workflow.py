"""Tests for entire workflows."""

from pathlib import Path

import matplotlib.colors
import matplotlib.testing.compare
import numpy
import pytest
import xarray as xr
from click.testing import CliRunner
from matplotlib import pyplot as plt

from wam2layers.cli import cli
from wam2layers.config import Config

try:  # Use cartopy if it's available
    import cartopy.feature as cfeature
    from cartopy import crs
except ImportError:
    crs = None
    cfeature = None


CFG_FILE_BACKWARD = Path("tests/test_data/config_rhine.yaml")
CFG_FILE_FORWARD = Path("tests/test_data/config_west_africa.yaml")


## First we define fixtures to help us write the tests more concisely
# Temporary directories that are consistent between tests:
@pytest.fixture(scope="session")
def temp_directory(tmp_path_factory) -> Path:
    return tmp_path_factory.mktemp("data")


@pytest.fixture(scope="session")
def output_dir(temp_directory: Path) -> Path:
    output_dir = temp_directory / "output"
    output_dir.mkdir()
    return output_dir


@pytest.fixture(scope="session")
def preprocessed_dir(temp_directory: Path) -> Path:
    return temp_directory / "preprocessed"


@pytest.fixture(scope="session")
def figures_directory(output_dir: Path) -> Path:
    fig_dir = output_dir / "figures"
    fig_dir.mkdir()
    return fig_dir


@pytest.fixture(scope="session")
def test_config_backward(
    temp_directory: Path, preprocessed_dir: Path, output_dir: Path
) -> str:
    """Return the path to the config for the preprocess, forward, and visualize tests.

    The processed and output data folders of this config are modified to use temporary
    directories.
    """
    cfg = Config.from_yaml(CFG_FILE_BACKWARD)
    cfg.preprocessed_data_folder = preprocessed_dir
    cfg.output_folder = output_dir
    tmp_config = temp_directory / "config_rhine.yaml"
    cfg.to_file(tmp_config)
    return str(tmp_config)


@pytest.fixture(scope="session")
def test_config_forward(
    temp_directory: Path, preprocessed_dir: Path, output_dir: Path
) -> str:
    """Return the path to the config for the forward tracking test(s).

    The output data folder of this config is modified to use a temporary directory.
    """
    cfg = Config.from_yaml(CFG_FILE_FORWARD)
    cfg.output_folder = output_dir
    tmp_config = temp_directory / "config_west_africa.yaml"
    cfg.to_file(tmp_config)
    return str(tmp_config)


# CLI runs. Session scope means that these fixtures will only run a single time
@pytest.fixture(scope="session")
def preprocess(test_config_backward: str) -> None:
    runner = CliRunner()
    runner.invoke(
        cli, ["preprocess", "era5", test_config_backward], catch_exceptions=False
    )


@pytest.fixture(scope="session")
def backtrack(test_config_backward: str, preprocess) -> None:
    runner = CliRunner()
    runner.invoke(cli, ["backtrack", test_config_backward], catch_exceptions=False)


@pytest.fixture(scope="session")
def visualize(test_config_backward: str, backtrack) -> None:
    runner = CliRunner()
    runner.invoke(
        cli, ["visualize", "output", test_config_backward], catch_exceptions=False
    )


@pytest.fixture(scope="session")
def forwardtrack(test_config_forward: str) -> None:
    runner = CliRunner()
    runner.invoke(cli, ["forwardtrack", test_config_forward], catch_exceptions=False)


def plot_difference(
    expected_output: xr.DataArray, output: xr.DataArray, fpath: Path
) -> None:
    """Plot the difference between the expected output and actual output."""
    if not fpath.parent.exists():
        fpath.parent.mkdir()  # Make figures directory if it does not exist yet

    if crs is not None:
        subplot_kw = dict(projection=crs.PlateCarree())
    else:
        subplot_kw = None

    fig, ax = plt.subplots(
        figsize=(16, 10),
        subplot_kw=subplot_kw,
    )
    (expected_output - output).plot.pcolormesh(
        ax=ax,
        x="longitude",
        cmap="RdBu",
        norm=matplotlib.colors.TwoSlopeNorm(vcenter=0),  # normalize cmap so 0 is center
    )
    if cfeature is not None:
        ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
        ax.add_feature(cfeature.BORDERS, linestyle="-", linewidth=0.2)

    ax.set_title(f"Difference between expected and actual '{output.name}' values.")
    fig.savefig(fpath, dpi=200)


def different_output_message(fig_path, output_path, expected_path):
    return (
        "Output results are different! Please verify your results. \n"
        f"A difference plot is available under {fig_path}. \n"
        "If you want to keep the new results and include it in your commit, \n"
        "you can update the reference data with the following command: \n"
        f"`cp {output_path} {expected_path}`"
    )


# Tests start now
class TestPreprocess:
    @pytest.fixture(autouse=True)
    def run_preprocess(self, preprocess):
        pass

    def test_preprocess(self, preprocessed_dir: Path, figures_directory: Path):
        output_path = preprocessed_dir / "2022-08-31_fluxes_storages.nc"
        assert output_path.exists()
        # verify outputs
        expected_path = Path(
            "tests/test_data/verify_output/2022-08-31_fluxes_storages.nc"
        )
        expected_output = xr.open_dataset(expected_path)
        output = xr.open_dataset(output_path)

        diff_figure = figures_directory / "preprocess_precip_diff.png"
        plot_difference(
            expected_output["precip"].sum("time"),
            output["precip"].sum("time"),
            diff_figure,
        )

        numpy.testing.assert_allclose(
            expected_output["precip"].values,
            output["precip"].values,
            err_msg=different_output_message(diff_figure, output_path, expected_path),
        )

    def test_log_file(self, preprocessed_dir):
        log_files = [i for i in preprocessed_dir.glob("wam2layers_*.log")]
        assert len(log_files) >= 1

    def test_config(self, preprocessed_dir):
        config_path = preprocessed_dir / "config_rhine.yaml"
        assert config_path.exists()


class TestBacktrack:
    @pytest.fixture(autouse=True)
    def run_backtrack(self, backtrack):
        pass

    def test_backtrack(self, output_dir: Path, figures_directory: Path):
        output_path = output_dir / "backtrack_2022-08-31T18-00.nc"
        assert output_path.exists()

        # verify outputs
        expected_path = Path(
            "tests/test_data/verify_output/backtrack_2022-08-31T18-00.nc"
        )
        expected_output = xr.open_dataset(expected_path)
        output = xr.open_dataset(output_path)

        diff_figure = figures_directory / "backtrack_diff.png"
        plot_difference(expected_output["e_track"], output["e_track"], diff_figure)

        numpy.testing.assert_allclose(
            expected_output["e_track"].values,
            output["e_track"].values,
            err_msg=different_output_message(diff_figure, output_path, expected_path),
        )

    def test_log_file(self, output_dir):
        log_files = [i for i in output_dir.glob("wam2layers_*.log")]
        len(log_files) >= 1

    def test_config(self, output_dir):
        config_path = output_dir / "config_rhine.yaml"
        assert config_path.exists()


class TestForwardtrack:
    @pytest.fixture(autouse=True)
    def run_forwardtrack(self, forwardtrack):
        pass

    def test_forwardtrack(self, output_dir: Path, figures_directory: Path):
        output_path = output_dir / "forwardtrack_1979-07-01T23-00.nc"
        assert output_path.exists()
        # verify outputs
        expected_path = Path(
            "tests/test_data/verify_output/forwardtrack_1979-07-01T23-00.nc"
        )
        expected_output = xr.open_dataset(expected_path)
        output = xr.open_dataset(output_path)

        for var in ["p_track_upper", "p_track_lower"]:
            diff_figure = figures_directory / f"forwardtrack_diff_{var}.png"
            expected_output[var] = expected_output[var]
            plot_difference(expected_output[var], output[var], diff_figure)

            numpy.testing.assert_allclose(
                expected_output[var].values,
                output[var].values,
                err_msg=different_output_message(
                    diff_figure, output_path, expected_path
                ),
            )

    def test_log_file(self, output_dir):
        log_files = [i for i in output_dir.glob("wam2layers_*.log")]
        assert len(log_files) >= 1

    def test_config(self, output_dir):
        config_path = output_dir / "config_west_africa.yaml"
        assert config_path.exists()


class TestVisualize:
    @pytest.fixture(autouse=True)
    def run_visualize(self, visualize):
        pass

    def test_visualize(self, output_dir: Path, figures_directory: Path):
        output_path = figures_directory / "cumulative_sources.png"
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
                f"`cp {output_path} {expected_output}`"
            )

    def test_log_file(self, output_dir):
        log_files = [i for i in output_dir.glob("wam2layers_*.log")]
        assert len(log_files) >= 2

    def test_config(self, output_dir):
        config_path = output_dir / "config_rhine.yaml"
        assert config_path.exists()
