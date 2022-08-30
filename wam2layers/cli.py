"""Command-line interface for the WAM2layers model."""
import click
from wam2layers.tracking.backtrack import cli as backtrack_cli

# Visualization requires optional dependencies that may not be installed.
try:
    from wam2layers.analysis.visualization import cli as visualize_cli
    include_visualize = True
except ImportError:
    include_visualize = False
    click.echo(
        "To enable visualization options, install matplotlib and cartopy. "
    )


@click.group()
def cli():
    pass


cli.add_command(backtrack_cli, name="backtrack")

# Visualization requires optional dependencies that may not be installed.
if include_visualize:
    cli.add_command(visualize_cli, name="visualize")
