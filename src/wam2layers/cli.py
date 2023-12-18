"""Command-line interface for the WAM2layers model."""
import logging
from datetime import datetime
from pathlib import Path

import click

from wam2layers.analysis.visualization import cli as visualize_cli
from wam2layers.config import Config
from wam2layers.preprocessing.cli import cli as preproc_cli
from wam2layers.tracking.backtrack import cli as backtrack_cli

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


@click.group()
@click.option("-c", "--config", required=True)
@click.pass_context
def cli(ctx, config):
    """Command line interface to WAM2layers."""
    ctx.obj = config
    config = Config.from_yaml(config)
    # logger.warn(f"Output path is {config.output_folder}.")
    # logger.warn("xxxxxxxxxxxxxxxxxxxxx")
    setup_logging(config.output_folder)


cli.add_command(preproc_cli, name="preprocess")
cli.add_command(backtrack_cli, name="backtrack")
cli.add_command(visualize_cli, name="visualize")
