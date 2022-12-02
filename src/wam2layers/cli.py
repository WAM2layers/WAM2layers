"""Command-line interface for the WAM2layers model."""
import click

from wam2layers.analysis.visualization import cli as visualize_cli
from wam2layers.preprocessing.cli import cli as preproc_cli
from wam2layers.tracking.backtrack import cli as backtrack_cli


@click.group()
def cli():
    """Command line interface to WAM2layers."""
    pass


cli.add_command(preproc_cli, name="preprocess")
cli.add_command(backtrack_cli, name="backtrack")
cli.add_command(visualize_cli, name="visualize")
