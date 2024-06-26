"""Command-line interface for the WAM2layers model.

This module contains all the functionality that makes it possible to run wam2layers on the command line, like so:

    wam2layers preprocess era5 floodcase.yaml
    wam2layers track floodcase.yaml

et cetera. It is built with [click](https://click.palletsprojects.com/en/8.1.x/).
"""

import logging
import shutil
import sys
from datetime import datetime
from pathlib import Path
from typing import Union

import click

from wam2layers import __version__
from wam2layers.analysis import visualization
from wam2layers.config import Config
from wam2layers.download import download_from_doi
from wam2layers.example_cases import AVAILABLE_CASES
from wam2layers.preprocessing.shared import prep_experiment
from wam2layers.tracking.backtrack import run_experiment as run_backtrack_experiment
from wam2layers.tracking.forwardtrack import (
    run_experiment as run_forwardtrack_experiment,
)

logger = logging.getLogger(__name__)


def setup_logging(log_path: Union[str, Path], debug: bool):
    """Configure logging behaviour.

    Messages with logging.INFO (and higher) are written to stdout
    Messages with logging.DEBUG (and higher) are written to wam2layers.log

    Args:
        log_path: Folder in which the logging file should be created.
        debug: If DEBUG level messages should also be printed to the terminal.
    """
    # https://stackoverflow.com/a/58828499
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")

    file_handler = logging.FileHandler(
        Path(log_path, f"wam2layers_{timestamp}.log"), mode="w"
    )
    stream_handler = logging.StreamHandler()

    if debug:
        stream_handler.setLevel(logging.DEBUG)
    else:
        stream_handler.setLevel(logging.INFO)

    logging.basicConfig(
        level=logging.DEBUG,
        format="%(message)s",
        handlers=[file_handler, stream_handler],
    )


def _copy_config_yaml(yaml_path, target_path):
    """Copy config file to the output path."""
    shutil.copy(yaml_path, target_path)


def get_wam_version():
    """Get the version of WAM2layers"""
    pkg_dir = Path(__file__).parent.absolute()
    return (
        f"wam2layers {__version__} from {pkg_dir} "
        f"(python {sys.version_info.major}.{sys.version_info.minor})"
    )


@click.group()
@click.pass_context
@click.version_option(message=get_wam_version())
@click.option("--debug", is_flag=True, help="Print debug statements to terminal.")
def cli(ctx, debug):
    """Command line interface to WAM2layers."""
    ctx.obj = dict()
    ctx.obj["debug"] = debug


# Command line setup for tracking
@cli.command()
@click.pass_context
@click.argument("config_file", type=click.Path(exists=True))
def track(ctx, config_file):
    """Run WAM2layers tracking experiment.

    CONFIG_FILE: Path to WAM2layers experiment configuration file.

    Usage examples:

        \b
        - wam2layers track path/to/cases/era5_2021.yaml
    """
    config = Config.from_yaml(config_file)
    log_path = config.output_folder
    setup_logging(log_path, ctx.obj["debug"])
    _copy_config_yaml(config_file, log_path)
    logger.info("Welcome to WAM2layers.")
    logger.info(f"Starting {config.tracking_direction} tracking experiment.")
    if config.tracking_direction == "backward":
        run_backtrack_experiment(config_file)
    else:
        run_forwardtrack_experiment(config_file)


@cli.command(hidden=True)
def backtrack():
    msg = (
        "The `backtrack` and `forwardtrack` commands have been removed in favor of "
        "`track`.\nPlease specify the tracking direction in your config file and use "
        "the `track` command instead."
        "\n"
        "\nFor example, in your config.yaml file add:"
        "\n    # Tracking"
        "\n    tracking_direction: backward"
    )
    raise ValueError(msg)


@cli.command(hidden=True)
def forwardtrack():
    msg = (
        "The `backtrack` and `forwardtrack` commands have been removed in favor of "
        "`track`.\nPlease specify the tracking direction in your config file and use "
        "the `track` command instead."
        "\n"
        "\nFor example, in your config.yaml file add:"
        "\n    # Tracking"
        "\n    tracking_direction: forward"
    )
    raise ValueError(msg)


# Command line setup for preprocess
@click.group()
def preproc_cli():
    """Pre-process raw input data for tracking with WAM2layers"""
    pass


@preproc_cli.command()
@click.pass_context
@click.argument("config_file", type=click.Path(exists=True))
def era5(ctx, config_file):
    """Preprocess ERA5 data for WAM2layers tracking experiments.

    CONFIG_FILE: Path to WAM2layers experiment configuration file.

    Usage examples:

        \b
        - wam2layers preprocess era5 path/to/cases/era5_2021.yaml
    """
    log_path = Config.from_yaml(config_file).preprocessed_data_folder
    setup_logging(log_path, ctx.obj["debug"])
    _copy_config_yaml(config_file, log_path)
    logger.info("Welcome to WAM2layers.")
    logger.info("Starting preprocessing ERA5 data.")
    prep_experiment(config_file, data_source="ERA5")


# Command line setup for visualization
@click.group()
def visualize_cli():
    """Visualize input or output data of a WAM2layers experiment"""
    ...


@visualize_cli.command()
@click.pass_context
@click.argument("config_file", type=click.Path(exists=True))
def input(ctx, config_file):
    """Visualize input data for experiment.

    CONFIG_FILE: Path to WAM2layers experiment configuration file.

    Usage examples:

        \b
        - wam2layers visualize input path/to/cases/era5_2021.yaml
    """
    log_path = Config.from_yaml(config_file).output_folder
    setup_logging(log_path, ctx.obj["debug"])
    _copy_config_yaml(config_file, log_path)
    logger.info("Welcome to WAM2layers.")
    logger.info("Starting visualizing input data.")
    visualization.visualize_input_data(config_file)


@visualize_cli.command()
@click.pass_context
@click.argument("config_file", type=click.Path(exists=True))
def output(ctx, config_file):
    """Visualize output data for experiment.

    CONFIG_FILE: Path to WAM2layers experiment configuration file.

    Usage examples:

        \b
        - wam2layers visualize output path/to/cases/era5_2021.yaml
    """
    log_path = Config.from_yaml(config_file).output_folder
    setup_logging(log_path, ctx.obj["debug"])
    _copy_config_yaml(config_file, log_path)
    logger.info("Welcome to WAM2layers.")
    logger.info("Starting visualizing output data.")
    visualization.visualize_output_data(config_file)


@visualize_cli.command()
@click.pass_context
@click.argument("config_file", type=click.Path(exists=True))
def both(ctx, config_file):
    """Visualize both input and output data for experiment."""
    log_path = Config.from_yaml(config_file).output_folder
    setup_logging(log_path, ctx.obj["debug"])
    _copy_config_yaml(config_file, log_path)
    logger.info("Welcome to WAM2layers.")
    logger.info("Starting visualizing both input and output data.")
    visualization.visualize_both(config_file)


cli.add_command(preproc_cli, name="preprocess")
cli.add_command(visualize_cli, name="visualize")


# Command line setup for downloading input data
@cli.command()
@click.pass_context
@click.argument("case", type=click.Choice(AVAILABLE_CASES.keys()))
def download(ctx, case):
    """Download input data for (example) cases."""
    logging.basicConfig(level=logging.INFO)
    doi = AVAILABLE_CASES.get(case, None)

    if doi is None:
        raise ValueError(
            f"Cannot download input data for {case}. Choose from {AVAILABLE_CASES.keys()}"
        )

    download_from_doi(doi, name=case)


if __name__ == "__main__":
    cli()
