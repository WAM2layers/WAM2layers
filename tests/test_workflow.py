"""Tests for entire workflows."""
from pathlib import Path

import matplotlib.colors
import matplotlib.testing.compare
import numpy
import xarray as xr
from click.testing import CliRunner
from matplotlib import pyplot as plt

from wam2layers.cli import cli

try:  # Use cartopy if it's available
    import cartopy.feature as cfeature
    from cartopy import crs
except ImportError:
    crs = None
    cfeature = None


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


def test_preprocess():
    runner = CliRunner()
    runner.invoke(
        cli,
        ["preprocess", "era5", "tests/test_data/test_config.yaml"],
        catch_exceptions=False,
    )

    output_path = Path("tests/tmp/preprocessed_data/2022-08-31_fluxes_storages.nc")
    assert output_path.exists()
    # verify outputs
    expected_output = xr.open_dataset(
        "tests/test_data/verify_output/2022-08-31_fluxes_storages.nc"
    )
    output = xr.open_dataset(output_path)

    diff_figure_path = Path("tests/tmp/output_data/figures/preprocess_precip_diff.png")
    plot_difference(
        expected_output["precip"].sum("time"),
        output["precip"].sum("time"),
        diff_figure_path,
    )

    numpy.testing.assert_allclose(
        expected_output["precip"].values,
        output["precip"].values,
        err_msg=(
            "Output results are different! Please verify your results. \n"
            f"A difference plot is available under {diff_figure_path}. \n"
            "If you want to keep the new results and include it in your commit, \n"
            "you can update the reference data with the following command: \n"
            "`cp tests/tmp/preprocessed_data/2022-08-31_fluxes_storages.nc tests/test_data/verify_output/2022-08-31_fluxes_storages.nc`"
        ),
    )
    # check log file
    log_files = [
        i for i in Path("tests/tmp/preprocessed_data/").glob("wam2layers_*.log")
    ]

    # check config yaml
    config_path = Path("tests/tmp/preprocessed_data/test_config.yaml")
    assert config_path.exists()


def test_backtrack():
    runner = CliRunner()
    runner.invoke(
        cli, ["backtrack", "tests/test_data/test_config.yaml"], catch_exceptions=False
    )

    output_path = Path("tests/tmp/output_data/backtrack_2022-08-31T18-00.nc")
    assert output_path.exists()
    # verify outputs
    expected_output = xr.open_dataset(
        "tests/test_data/verify_output/backtrack_2022-08-31T18-00.nc"
    )
    output = xr.open_dataset(output_path)

    diff_figure_path = Path("tests/tmp/output_data/figures/backtrack_diff.png")
    plot_difference(expected_output["e_track"], output["e_track"], diff_figure_path)

    numpy.testing.assert_allclose(
        expected_output["e_track"].values,
        output["e_track"].values,
        err_msg=(
            "Output results are different! Please verify your results. \n"
            f"A difference plot is available under {diff_figure_path}. \n"
            "If you want to keep the new results and include it in your commit, \n"
            "you can update the reference data with the following command: \n"
            "`cp tests/tmp/output_data/backtrack_2022-08-31T18-00.nc tests/test_data/verify_output/backtrack_2022-08-31T18-00.nc`"
        ),
    )

    # check log file
    log_files = [i for i in Path("tests/tmp/output_data/").glob("wam2layers_*.log")]

    # check config yaml
    config_path = Path("tests/tmp/output_data/test_config.yaml")
    assert config_path.exists()


def test_forwardtrack():
    runner = CliRunner()
    runner.invoke(
        cli,
        ["forwardtrack", "tests/test_data/config_west_africa.yaml"],
        catch_exceptions=False,
    )

    output_path = Path("tests/tmp/output_data/forwardtrack_1979-07-01T23-00.nc")
    assert output_path.exists()
    # verify outputs
    expected_path = Path(
        "tests/test_data/verify_output/forwardtrack_1979-07-01T23-00.nc"
    )
    expected_output = xr.open_dataset(expected_path)
    output = xr.open_dataset(output_path)

    for var in ["p_track_upper", "p_track_lower"]:
        diff_figure_path = Path(
            f"tests/tmp/output_data/figures/forwardtrack_diff_{var}.png"
        )
        expected_output[var] = expected_output[var]
        plot_difference(expected_output[var], output[var], diff_figure_path)

        numpy.testing.assert_allclose(
            expected_output[var].values,
            output[var].values,
            err_msg=(
                "Output results are different! Please verify your results. \n"
                f"A difference plot is available under {diff_figure_path}. \n"
                "If you want to keep the new results and include it in your commit, \n"
                "you can update the reference data with the following command: \n"
                f"`cp {output_path} {expected_path}`"
            ),
        )

    # check log file
    log_files = [i for i in Path("tests/tmp/output_data/").glob("wam2layers_*.log")]

    # check config yaml
    config_path = Path("tests/tmp/output_data/config_west_africa.yaml")
    assert config_path.exists()


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

    # check log file
    log_files = [i for i in Path("tests/tmp/output_data/").glob("wam2layers_*.log")]

    # check config yaml
    config_path = Path("tests/tmp/output_data/test_config.yaml")
    assert config_path.exists()
