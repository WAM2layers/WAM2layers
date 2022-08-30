from cartopy import crs, feature as cfeature
import click
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from pathlib import Path
from wam2layers.preprocessing.preprocessing import get_grid_info
from wam2layers.tracking.backtrack import parse_config, input_path, output_path, load_region

def polish(ax, region):
    ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax.add_feature(cfeature.BORDERS, linestyle="-", linewidth=0.2)
    ax.set_xlim(region.longitude.min(), region.longitude.max())
    ax.set_ylim(region.latitude.min(), region.latitude.max())


def visualize_input_data(config_file):
    raise NotImplementedError()
    config = parse_config(config_file)
    for date in config["datelist"]:
        ds = xr.open_dataset(input_path(date, config))
        import IPython; IPython.embed(); quit()
    # do stuff


def visualize_output_data(config_file):
    raise NotImplementedError()
    config = parse_config(config_file)
    for date in config["datelist"]:
        ds = xr.open_dataset(output_path(date, config))
        import IPython; IPython.embed(); quit()
    # do stuff



def visualize_both(config_file):
    """Diagnostic figure with four subplots combining input and output data."""
    config = parse_config(config_file)
    region = load_region(config)
    a_gridcell, lx, ly = get_grid_info(region)

    out_dir = Path(config["output_folder"]) / "figures"
    out_dir.mkdir(exist_ok=True, parents=True)

    for date in config["datelist"]:
        print(date)
        ds_in = xr.open_dataset(input_path(date, config))
        ds_out = xr.open_dataset(output_path(date, config))

        # Make figure
        fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(2, 2, subplot_kw=dict(projection=crs.PlateCarree()), figsize=(14,8))

        ax1.set_title("Tracked precipitation" + date.strftime("%Y%m%d"))
        precip = ds_in.precip.sum('time') * region / a_gridcell[:, None] * 1000
        precip.plot.contourf(ax=ax1, cmap='Blues', levels=np.arange(0, 50, 5))
        polish(ax1)

        ax2.set_title("Moisture source")
        e_track = ds_out.e_track / a_gridcell[:, None] * 1000
        e_track.plot.contourf(ax=ax2, cmap='GnBu', levels=np.arange(0.0, 1.0, 0.1))
        polish(ax2)

        ax3.set_title("S track upper layer")
        s_track_upper = ds_out.s_track_upper / a_gridcell[:, None] * 1000
        s_track_upper.plot.contourf(ax=ax3, cmap="YlOrRd", levels=np.arange(0.0, 5, 0.5))
        ds_in.mean('time').plot.streamplot(x="longitude", y="latitude", u="fx_upper", v="fy_upper", ax=ax3, color="black")
        polish(ax3)

        ax4.set_title("S track lower layer")
        s_track_lower = ds_out.s_track_lower / a_gridcell[:, None] * 1000
        s_track_lower.plot.contourf(ax=ax4, cmap="YlOrRd", levels=np.arange(0.0, 5, 0.5))
        ds_in.mean('time').plot.streamplot(x="longitude", y="latitude", u="fx_lower", v="fy_lower", ax=ax4, color="black")
        polish(ax4)

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
    pass


@cli.command()
@click.argument('config_file', type=click.Path(exists=True))
def input(config_file):
    """Visualize input data for experiment."""
    visualize_input_data(config_file)


@cli.command()
@click.argument('config_file', type=click.Path(exists=True))
def output(config_file):
    """Visualize output data for experiment."""
    visualize_output_data(config_file)


@cli.command()
@click.argument('config_file', type=click.Path(exists=True))
def both(config_file):
    """Visualize both input and output data for experiment."""
    visualize_both(config_file)

if __name__ == "__main__":
    cli()
