import logging
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
from cmocean import cm

from wam2layers.config import Config
from wam2layers.tracking.io import input_path, load_tagging_region, output_path

logger = logging.getLogger(__name__)


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


def _plot_input(config: Config, ax):
    """
    Return subplot with input.
    backtrack = precipitation
    forwardtrack = evaporation
    """
    # Load config and some useful stuf.
    region = load_tagging_region(config)

    # Load data
    start = config.tagging_start_date
    end = config.tagging_end_date
    dates = pd.date_range(
        start=start,
        end=end,
        freq=config.output_frequency,
        inclusive="left",
    )

    input_files = []
    for date in dates[:-1]:
        input_files.append(input_path(date, config))
    ds = xr.open_mfdataset(input_files, combine="nested", concat_dim="time")
    if config.tracking_direction == "backward":
        subset = ds.precip.sel(time=slice(start, end))
    else:
        subset = ds.evap.sel(time=slice(start, end))

    # TODO: check if backtrack files or forwardtrack files were in output
    input = (subset * region * 3600).sum("time").compute()

    # Make figure
    input.plot(ax=ax, cmap=cm.rain, cbar_kwargs=dict(fraction=0.05, shrink=0.5))
    ax.set_title("Cumulative input during tagging [mm]", loc="left")
    polish(ax, region.where(region > 0, drop=True))


def _plot_output(config: Config, ax):
    """
    Return subplot with tracked moisture
    forwardtrack = precipitation
    backtrack = evaporation
    """
    region = load_tagging_region(config)

    # Load data
    dates = pd.date_range(
        start=config.tracking_start_date,
        end=config.tracking_end_date,
        freq=config.output_frequency,
        inclusive="left",
    )

    if config.tracking_direction == "backward":
        output_files = [output_path(date, config) for date in dates[:-1]]
    else:
        output_files = [output_path(date, config) for date in dates[1:]]

    if not all([Path(file).exists() for file in output_files]):
        raise FileNotFoundError(f"Could not find all files: {output_files}.")
    elif len(output_files) == 0:
        raise ValueError("Too few output files to visualize.")

    ds = xr.open_mfdataset(output_files, combine="nested", concat_dim="time")
    if config.tracking_direction == "backward":
        out_track = ds.e_track.sum("time")
    else:
        out_track = ds.p_track_lower.sum("time") + ds.p_track_upper.sum("time")

    # Make figure
    out_track.plot(
        ax=ax,
        vmin=0,
        robust=True,
        cmap=cm.rain,
        cbar_kwargs=dict(fraction=0.05, shrink=0.5),
    )

    ax.set_title("Accumulated tracked moisture [mm]", loc="left")

    # Add source region outline
    # TODO: fix this, because no region contours are shown
    region.plot.contour(ax=ax, levels=[1], colors="k")
    polish(ax, region)


def visualize_input_data(config_file):
    """An figure showing the cumulative moisture inputs.

    TODO: make figure creation independent of case.
    """
    config = Config.from_yaml(config_file)

    # Make figure
    fig = plt.figure(figsize=(16, 10))
    ax = fig.add_subplot(111, projection=crs.PlateCarree())

    _plot_input(config, ax)

    # Save
    out_dir = Path(config.output_folder) / "figures"
    out_dir.mkdir(exist_ok=True, parents=True)
    path_to_fig = out_dir / "input_event.png"
    fig.savefig(path_to_fig, dpi=200)
    logger.info(f"Figure of cumulative moisture inputs written to {path_to_fig}.")


def visualize_output_data(config_file):
    """An figure showing the cumulative moisture origins.

    TODO: make figure creation independent of case.
    """
    # Load config and some useful stuf.
    config = Config.from_yaml(config_file)

    # Make figure
    fig = plt.figure(figsize=(16, 10))
    ax = fig.add_subplot(111, projection=crs.PlateCarree())
    _plot_output(config, ax)

    # Save
    out_dir = Path(config.output_folder) / "figures"
    out_dir.mkdir(exist_ok=True, parents=True)
    path_to_fig = out_dir / "cumulative_sources.png"
    fig.savefig(path_to_fig, dpi=200)
    logger.info(f"Figure of cumulative moisture origins written to {path_to_fig}.")


def visualize_both(config_file):
    """Diagnostic figure with four subplots combining input and output data."""
    # Load config and some useful stuf.
    config = Config.from_yaml(config_file)

    # Make figure
    fig, [ax1, ax2] = plt.subplots(
        2, 1, figsize=(16, 10), subplot_kw=dict(projection=crs.PlateCarree())
    )

    _plot_input(config, ax1)
    _plot_output(config, ax2)

    # Save
    out_dir = Path(config.output_folder) / "figures"
    out_dir.mkdir(exist_ok=True, parents=True)
    path_to_fig = out_dir / "summary_subplots.png"
    fig.savefig(path_to_fig, dpi=200)
    logger.info(f"Diagnostic figure of input and output data written to {path_to_fig}.")


# TODO: this part seems broken
def visualize_snapshots(config_file):
    """Diagnostic figure with four subplots combining input and output data."""
    # TODO: make this work for forwardtrack
    config = Config.from_yaml(config_file)
    dates = pd.date_range(
        start=config.tracking_start_date,
        end=config.tracking_end_date,
        freq=config.output_frequency,
        inclusive="left",
    )
    region = load_tagging_region(config)

    out_dir = Path(config.output_folder) / "figures"
    out_dir.mkdir(exist_ok=True, parents=True)

    for date in dates:
        logger.info(date)
        ds_in = xr.open_dataset(input_path(date, config))
        ds_out = xr.open_dataset(output_path(date, config))

        # Make figure
        fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(
            2, 2, subplot_kw=dict(projection=crs.PlateCarree()), figsize=(14, 8)
        )

        ax1.set_title("Tagged moisture" + date.strftime("%Y%m%d"))
        precip = ds_in.precip.sum("time") * region
        precip.plot(ax=ax1, cmap="Blues")
        polish(ax1, region)

        ax2.set_title("Tracked moisture")
        e_track = ds_out.e_track
        e_track.plot(ax=ax2, cmap="GnBu")
        polish(ax2, region)

        ax3.set_title("S track upper layer")
        s_track_upper = ds_out.s_track_upper_restart
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
        s_track_lower = ds_out.s_track_lower_restart
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
        logger.info(f"Snapshot figure written to {output_file}.")
