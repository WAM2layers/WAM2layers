from pathlib import Path

import click
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
from cmocean import cm

from wam2layers.preprocessing.shared import get_grid_info
from wam2layers.tracking.backtrack import (input_path, load_region,
                                           output_path, parse_config)


def try_import_cartopy():
    """Import cartopy if it is available; else raise."""
    from importlib import import_module

    global crs
    global cfeature

    try:
        crs = import_module("cartopy.crs")
        cfeature = import_module("cartopy.feature")
    except ImportError as exec:
        message = "This function requires cartopy. Cartopy is most easily installed with conda. Please refer to the documentation."
        raise ImportError(message) from exec


def polish(ax, region):
    ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax.add_feature(cfeature.BORDERS, linestyle="-", linewidth=0.2)
    ax.set_xlim(region.longitude.min(), region.longitude.max())
    ax.set_ylim(region.latitude.min(), region.latitude.max())


def _plot_precip(config, ax):
    """Return subplot with precip."""
    # Load config and some usful stuf.
    region = load_region(config)

    # Load data
    dates = pd.date_range(
        start=config["preprocess_start_date"],
        end=config["preprocess_end_date"],
        freq=config[
            "output_frequency"
        ],  # Should be output frequency, since this is used to save the data
        inclusive="left",
    )

    input_files = []
    for date in dates:
        input_files.append(input_path(date, config))

    ds = xr.open_mfdataset(input_files, combine="nested", concat_dim="time")
    # TODO: make region time-dependent
    start = config["event_start_date"]
    end = config["event_end_date"]
    subset = ds.precip.sel(time=slice(start, end))
    precip = (subset * region * 3600).sum("time").compute()

    # Make figure
    precip.plot(ax=ax, cmap=cm.rain, cbar_kwargs=dict(fraction=0.05, shrink=0.5))
    ax.set_title("Cumulative precipitation during event [mm]", loc="left")
    polish(ax, region.where(region > 0, drop=True))


def _plot_evap(config, ax):
    """Return subplot with tracked evaporation."""
    region = load_region(config)
    a_gridcell, lx, ly = get_grid_info(region)

    # Load data
    dates = pd.date_range(
        start=config["track_start_date"],
        end=config["track_end_date"],
        freq=config["output_frequency"],
        inclusive="left",
    )

    output_files = []
    for date in dates:
        output_files.append(output_path(date, config))

    ds = xr.open_mfdataset(output_files, combine="nested", concat_dim="time")
    e_track = ds.e_track.sum("time").compute() * 1000 / a_gridcell[:, None]

    # Make figure
    e_track.plot(
        ax=ax, vmin=0, cmap=cm.rain, cbar_kwargs=dict(fraction=0.05, shrink=0.5)
    )
    e_track.plot.contour(ax=ax, levels=[0.1, 1], colors=["lightgrey", "grey"])
    ax.set_title("Moisture source of extreme precip [mm]", loc="left")

    # Add source region outline
    region.plot.contour(ax=ax, levels=[1], colors="k")
    polish(ax, region)


def visualize_input_data(config_file):
    """An figure showing the cumulative moisture inputs.

    TODO: make figure creation independent of case.
    """
    config = parse_config(config_file)

    # Make figure
    fig = plt.figure(figsize=(16, 10))
    ax = fig.add_subplot(111, projection=crs.PlateCarree())

    _plot_precip(config, ax)

    # Save
    out_dir = Path(config["output_folder"]) / "figures"
    out_dir.mkdir(exist_ok=True, parents=True)
    fig.savefig(out_dir / "input_event.png", dpi=200)


def visualize_output_data(config_file):
    """An figure showing the cumulative moisture origins.

    TODO: make figure creation independent of case.
    """
    # Load config and some usful stuf.
    config = parse_config(config_file)

    # Make figure
    fig = plt.figure(figsize=(16, 10))
    ax = fig.add_subplot(111, projection=crs.PlateCarree())

    _plot_evap(config, ax)

    # Save
    out_dir = Path(config["output_folder"]) / "figures"
    out_dir.mkdir(exist_ok=True, parents=True)
    fig.savefig(out_dir / "cumulative_sources.png", dpi=200)


def visualize_both(config_file):
    """Diagnostic figure with four subplots combining input and output data."""
    # Load config and some usful stuf.
    config = parse_config(config_file)

    # Make figure
    fig, [ax1, ax2] = plt.subplots(
        2, 1, figsize=(16, 10), subplot_kw=dict(projection=crs.PlateCarree())
    )

    _plot_precip(config, ax1)
    _plot_evap(config, ax2)

    # Save
    out_dir = Path(config["output_folder"]) / "figures"
    out_dir.mkdir(exist_ok=True, parents=True)
    fig.savefig(out_dir / "summary_subplots.png", dpi=200)


def visualize_snapshots(config_file):
    """Diagnostic figure with four subplots combining input and output data."""
    config = parse_config(config_file)
    dates = pd.date_range(
        start=config["track_start_date"],
        end=config["track_end_date"],
        freq=config["output_frequency"],
        inclusive="left",
    )
    region = load_region(config)
    a_gridcell, lx, ly = get_grid_info(region)

    out_dir = Path(config["output_folder"]) / "figures"
    out_dir.mkdir(exist_ok=True, parents=True)

    for date in dates:
        print(date)
        ds_in = xr.open_dataset(input_path(date, config))
        ds_out = xr.open_dataset(output_path(date, config))

        # Make figure
        fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(
            2, 2, subplot_kw=dict(projection=crs.PlateCarree()), figsize=(14, 8)
        )

        ax1.set_title("Tracked precipitation" + date.strftime("%Y%m%d"))
        precip = ds_in.precip.sum("time") * region / a_gridcell[:, None] * 1000
        precip.plot(ax=ax1, cmap="Blues")
        polish(ax1, region)

        ax2.set_title("Moisture source")
        e_track = ds_out.e_track / a_gridcell[:, None] * 1000
        e_track.plot(ax=ax2, cmap="GnBu")
        polish(ax2, region)

        ax3.set_title("S track upper layer")
        s_track_upper = ds_out.s_track_upper_restart / a_gridcell[:, None] * 1000
        s_track_upper.plot(ax=ax3, cmap="YlOrRd")
        ds_in.mean("time").plot.streamplot(
            x="longitude",
            y="latitude",
            u="fx_upper",
            v="fy_upper",
            ax=ax3,
            color="black",
        )
        polish(ax3, region)

        ax4.set_title("S track lower layer")
        s_track_lower = ds_out.s_track_lower_restart / a_gridcell[:, None] * 1000
        s_track_lower.plot(ax=ax4, cmap="YlOrRd")
        ds_in.mean("time").plot.streamplot(
            x="longitude",
            y="latitude",
            u="fx_lower",
            v="fy_lower",
            ax=ax4,
            color="black",
        )
        polish(ax4, region)

        # Save
        output_file = out_dir / f"input_output_{date.strftime('%Y%m%d')}.png"
        fig.savefig(output_file)
        plt.close()


###########################################################################
# The code below makes it possible to run wam2layers from the command line:
# >>> python visualization.py input path/to/cases/era5_2021.yaml
# or even:
# >>> wam2layers visualize output path/to/cases/era5_2021.yaml
###########################################################################


@click.group()
def cli():
    """Visualize input or output data of a WAM2layers experiment"""
    pass


@cli.command()
@click.argument("config_file", type=click.Path(exists=True))
def input(config_file):
    """Visualize input data for experiment."""
    try_import_cartopy()
    visualize_input_data(config_file)


@cli.command()
@click.argument("config_file", type=click.Path(exists=True))
def output(config_file):
    """Visualize output data for experiment."""
    try_import_cartopy()
    visualize_output_data(config_file)


@cli.command()
@click.argument("config_file", type=click.Path(exists=True))
def both(config_file):
    """Visualize both input and output data for experiment."""
    try_import_cartopy()
    visualize_both(config_file)


@cli.command()
@click.argument("config_file", type=click.Path(exists=True))
def snapshots(config_file):
    """Visualize input and output snapshots for experiment."""
    try_import_cartopy()
    visualize_snapshots(config_file)


if __name__ == "__main__":
    cli()
