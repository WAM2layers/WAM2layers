"""Command-line interface for the WAM2layers model.

This module contains all the functionality that makes it possible to run wam2layers on the command line, like so:

    wam2layers preprocess era5 floodcase.yaml
    wam2layers track floodcase.yaml

et cetera. It is built with [click](https://click.palletsprojects.com/en/8.1.x/).
"""
import logging
import shutil
from datetime import datetime
from pathlib import Path
import sys

import click

from wam2layers import __version__
from wam2layers.analysis import visualization
from wam2layers.config import Config
from wam2layers.preprocessing.era5 import prep_experiment
from wam2layers.tracking.backtrack import run_experiment as run_backtrack_experiment
from wam2layers.tracking.forwardtrack import (
    run_experiment as run_forwardtrack_experiment,
)

logger = logging.getLogger(__name__)


def setup_logging(log_path):
    """Configure logging behaviour.

    Messages with logger.INFO (and higher) are written to stdout
    Messages with logger.DEBUG (and higher) are written to wam2layers.log
    """
    # https://stackoverflow.com/a/58828499
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")

    file_handler = logging.FileHandler(
        Path(log_path, f"wam2layers_{timestamp}.log"), mode="w"
    )
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    logging.basicConfig(
        level=logging.INFO,
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
@click.version_option(message=get_wam_version())
def cli():
    """Command line interface to WAM2layers."""
    pass


# Command line setup for tracking


@cli.command()
@click.argument("config_file", type=click.Path(exists=True))
def track(config_file):
    """Run WAM2layers tracking experiment.

    CONFIG_FILE: Path to WAM2layers experiment configuration file.

    Usage examples:

        \b
        - wam2layers track path/to/cases/era5_2021.yaml
    """
    config = Config.from_yaml(config_file)
    log_path = config.output_folder
    setup_logging(log_path)
    _copy_config_yaml(config_file, log_path)
    logger.info("Welcome to WAM2layers.")
    logger.info(f"Starting {config.tracking_direction} tracking experiment.")
    if config.tracking_direction == "backward":
        run_backtrack_experiment(config_file)
    else:
        run_forwardtrack_experiment(config_file)


@cli.command()
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


@cli.command()
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
@click.argument("config_file", type=click.Path(exists=True))
def era5(config_file):
    """Preprocess ERA5 data for WAM2layers tracking experiments.

    CONFIG_FILE: Path to WAM2layers experiment configuration file.

    Usage examples:

        \b
        - wam2layers preprocess era5 path/to/cases/era5_2021.yaml
    """
    log_path = Config.from_yaml(config_file).preprocessed_data_folder
    setup_logging(log_path)
    _copy_config_yaml(config_file, log_path)
    logger.info("Welcome to WAM2layers.")
    logger.info("Starting preprocessing ERA5 data.")
    prep_experiment(config_file)


# Command line setup for visualization
@click.group()
def visualize_cli():
    """Visualize input or output data of a WAM2layers experiment"""
    ...


@visualize_cli.command()
@click.argument("config_file", type=click.Path(exists=True))
def input(config_file):
    """Visualize input data for experiment.

    CONFIG_FILE: Path to WAM2layers experiment configuration file.

    Usage examples:

        \b
        - wam2layers visualize input path/to/cases/era5_2021.yaml
    """
    log_path = Config.from_yaml(config_file).output_folder
    setup_logging(log_path)
    _copy_config_yaml(config_file, log_path)
    logger.info("Welcome to WAM2layers.")
    logger.info("Starting visualizing input data.")
    visualization.visualize_input_data(config_file)


@visualize_cli.command()
@click.argument("config_file", type=click.Path(exists=True))
def output(config_file):
    """Visualize output data for experiment.

    CONFIG_FILE: Path to WAM2layers experiment configuration file.

    Usage examples:

        \b
        - wam2layers visualize output path/to/cases/era5_2021.yaml
    """
    log_path = Config.from_yaml(config_file).output_folder
    setup_logging(log_path)
    _copy_config_yaml(config_file, log_path)
    logger.info("Welcome to WAM2layers.")
    logger.info("Starting visualizing output data.")
    visualization.visualize_output_data(config_file)


@visualize_cli.command()
@click.argument("config_file", type=click.Path(exists=True))
def both(config_file):
    """Visualize both input and output data for experiment."""
    log_path = Config.from_yaml(config_file).output_folder
    setup_logging(log_path)
    _copy_config_yaml(config_file, log_path)
    logger.info("Welcome to WAM2layers.")
    logger.info("Starting visualizing both input and output data.")
    visualization.visualize_both(config_file)


@visualize_cli.command()
@click.argument("config_file", type=click.Path(exists=True))
def snapshots(config_file):
    """Visualize input and output snapshots for experiment."""
    log_path = Config.from_yaml(config_file).output_folder
    setup_logging(log_path)
    _copy_config_yaml(config_file, log_path)
    logger.info("Welcome to WAM2layers.")
    logger.info("Starting visualizing input and output snapshots.")
    visualization.visualize_snapshots(config_file)


cli.add_command(preproc_cli, name="preprocess")
cli.add_command(visualize_cli, name="visualize")


if __name__ == "__main__":
    cli()
