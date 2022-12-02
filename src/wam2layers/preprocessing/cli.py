import click

from wam2layers.preprocessing.era5 import cli as era5_cli


@click.group()
def cli():
    """Pre-process raw input data for tracking with WAM2layers"""
    pass


cli.add_command(era5_cli, name="era5")
