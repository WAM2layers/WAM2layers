"""Command-line interface for the WAM2layers model."""
import logging
from datetime import datetime

import click

from wam2layers.analysis.visualization import cli as visualize_cli
from wam2layers.preprocessing.cli import cli as preproc_cli
from wam2layers.tracking.backtrack import cli as backtrack_cli


def setup_logging():
    """Configure logging behaviour.

    Messages with logger.INFO (and higher) are written to stdout
    Messages with logger.DEBUG (and higher) are written to wam2layers.log
    """
    # https://stackoverflow.com/a/58828499
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")

    file_handler = logging.FileHandler(f"wam2layers_{timestamp}.log", mode="w")
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(message)s",
        handlers=[file_handler, stream_handler],
    )


@click.group()
def cli():
    """Command line interface to WAM2layers."""
    setup_logging()


cli.add_command(preproc_cli, name="preprocess")
cli.add_command(backtrack_cli, name="backtrack")
cli.add_command(visualize_cli, name="visualize")
