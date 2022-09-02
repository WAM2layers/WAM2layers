from pathlib import Path

import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter

from preprocessing.preprocessing import get_grid_info


def make_diagnostic_figures(
    date,
    region,
    fx_upper,
    fy_upper,
    fx_lower,
    fy_lower,
    precip,
    s_track_upper_mean,
    s_track_lower_mean,
    e_track,
):
    """Visualize fields during the simulation."""
    # import IPython; IPython.embed(); quit()
    output_dir = Path("../figures")
    output_dir.mkdir(exist_ok=True)

    # Load data
    u = xr.open_dataset("/data/volume_2/era5_2021/FloodCase_202107_u.nc")

    # Get grid info
    lat = u.latitude.values
    lon = u.longitude.values
    a_gridcell, l_ew_gridcell, l_mid_gridcell = get_grid_info(u)

    # TODO: improve this
    a_gridcell = a_gridcell[:, None]

    my_projection = ccrs.PlateCarree(central_longitude=0)

    def polish(ax):
        ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.8)
        ax.add_feature(cartopy.feature.BORDERS, linestyle="-", linewidth=0.2)

        ax.set_xticks(np.arange(-180, 181, 20), crs=my_projection)
        ax.set_yticks(np.arange(-90, 91, 20), crs=my_projection)
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)
        ax.contour(lon, lat, region)
        ax.set_xlim(-50, 30)
        ax.set_ylim(30, 60)

    precip_track = np.arange(0.0, 50.0, 5)
    S_track = np.arange(0.0, 5, 0.5)
    E_track = np.arange(0.0, 1.0, 0.1)

    fig1 = plt.figure(figsize=(14, 8))

    ax1 = plt.subplot(221, projection=my_projection)
    cb1 = ax1.contourf(
        lon,
        lat,
        (precip.sum(axis=0) * region / a_gridcell) * 1000,
        precip_track,
        cmap=plt.cm.Blues,
        extend="max",
    )  # We plot a colormesh using the gist_ncar colormap.
    polish(ax1)
    ax1.set_title("Tracked precipitation" + date.strftime("%Y%m%d"))

    ax2 = plt.subplot(222, projection=my_projection)
    cb2 = ax2.contourf(
        lon, lat, (e_track / a_gridcell) * 1000, E_track, cmap=plt.cm.GnBu, extend="max"
    )  # We plot a colormesh using the gist_ncar colormap.
    polish(ax2)
    ax2.set_title("Moisture source")

    ax3 = plt.subplot(223, projection=my_projection)
    cb3 = ax3.contourf(
        lon,
        lat,
        (s_track_upper_mean / a_gridcell) * 1000,
        S_track,
        cmap=plt.cm.YlOrRd,
        extend="max",
    )  # We plot a colormesh using the gist_ncar colormap.
    ax3.quiver(
        lon[::5],
        lat[::5],
        fx_upper.mean(axis=0)[::5, ::5],
        fy_upper.mean(axis=0)[::5, ::5],
        color="black",
        width=0.003,
        alpha=0.5,
    )
    polish(ax3)
    ax3.set_title("S track upper layer")

    ax4 = plt.subplot(224, projection=my_projection)
    cb3 = ax4.contourf(
        lon,
        lat,
        (s_track_lower_mean / a_gridcell) * 1000,
        S_track,
        cmap=plt.cm.YlOrRd,
        extend="max",
    )  # We plot a colormesh using the gist_ncar colormap.
    ax4.quiver(
        lon[::5],
        lat[::5],
        fx_lower.mean(axis=0)[::5, ::5],
        fy_lower.mean(axis=0)[::5, ::5],
        color="black",
        width=0.003,
        alpha=0.5,
    )
    polish(ax4)
    ax4.set_title("S track lower layer")

    new_axis = fig1.add_axes([0.55, 0.50, 0.35, 0.015])  # left, bottom, width, height
    fig1.colorbar(cb2, cax=new_axis, orientation="horizontal")

    new_axis2 = fig1.add_axes([0.15, 0.50, 0.35, 0.015])
    fig1.colorbar(cb1, cax=new_axis2, orientation="horizontal")

    new_axis3 = fig1.add_axes([0.30, 0.08, 0.35, 0.015])
    fig1.colorbar(cb3, cax=new_axis3, orientation="horizontal")
    plt.savefig(output_dir / f"tracking_{date.strftime('%Y%m%d')}.png")
    plt.close()
