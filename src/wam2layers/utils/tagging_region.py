"""
WAM2layers Tagging Region Designer functions

@author: Vincent de Feiter | GitHub: vincentdefeiter

Adapted by Chris Weijenborg

"""
import logging

# Import
import cartopy.crs as ccrs  # used for plotting on a map
import geopandas as gpd  # reading shapefiles
import matplotlib.patheffects as pe  # additional plotting functionalities
import matplotlib.pyplot as plt  # used for plotting purposes
import numpy as np  # used for array calculations
import regionmask  # used for region mapping
import xarray as xr  # used to read .nc files
from cartopy import feature as cfeature

logger = logging.getLogger(__name__)


def mask_with_regionmask(region_source, regions, input_file, output_dir, return_plot):
    """tagging region by using a named region from regionmask.

    This function builds a tagging region based on available regions from the Region Mask Python Package.

    Args:
        region_source (str): indicate here a region source from the Region Mask Python Package (Giorgi, SREX, AR6, PRUDENCE).
        regions (list): single or multiple regions as list.
        reference_file: a reference file of the preprocessing which is used to
            extract the dimensions of the data used.
        output_dir: path where to store the created tagging region.
        return_plot: if true, return a plot that shows the tagging region on a map.

    Returns:
        axes on which the region is plotted if `return_plot` is True

    """

    # Retrieve file dimensions
    file = xr.open_dataset(input_file)

    longitude = np.array(file.longitude)
    latitude = np.array(file.latitude)

    # Set available regions
    if "ar6" in region_source:
        name, region_source = region_source.split(".")
        try:
            available_regions = getattr(regionmask.defined_regions.ar6, region_source)

        except:
            raise KeyError(
                'region_source not specified correctly, please check the "regionmask.defined_regions" options in the Region Mask documentation for more details'
            )
    else:
        try:
            available_regions = getattr(regionmask.defined_regions, region_source)
        except:
            raise KeyError(
                'region_source not specified correctly, please check the "regionmask.defined_regions" options in the Region Mask documentation for more details'
            )

    # create (a) mask(s) of the required region(s):
    mask = getattr(available_regions, "mask_3D")(longitude, latitude)

    # Specify region ID - what to store
    region_ids = np.array(regions)

    if (
        "srex" in region_source
        or "giorgi" in region_source
        or "prudence" in region_source
    ):
        region_ids = region_ids - 1
    if "ocean" in region_source or "land" in region_source:
        raise KeyError(
            "Entry for ar6 not recognized. Please only use the ar6.all options to make the selection. This option ensures proper indexing of the region."
        )

    region_ids = region_ids.tolist()

    if type(region_ids) == int:
        raise KeyError(
            "regions not specified correctly, please check whether your input is provided as a list [...,...], even if it is a single tagging region!"
        )

    new_masks = mask.isel(region=region_ids).values

    # Make WAM2layers fitting
    new_masks = np.where(
        new_masks == False, 0, 1
    )  # Replace True and False by 1 and 0's
    new_masks = np.nanmean(new_masks, axis=0)  # Multiple masks into 1D array

    # Export as tagging region for WAM2layers

    # create xarray dataset
    data = new_masks

    export = xr.Dataset(
        {"tagging_region": (["latitude", "longitude"], data.astype(float))}
    )

    # set coordinates
    export["longitude"] = ("longitude", longitude)
    export["latitude"] = ("latitude", latitude)

    # save to NetCDF file
    export.to_netcdf(output_dir)

    logger.info(f"Stored tagging region in {output_dir}.")

    # Plot as a check
    if return_plot:
        # Visualise all regions available

        fontsizes = 10
        pads = 20

        text_kws = dict(
            bbox=dict(color="none"),
            path_effects=[pe.withStroke(linewidth=2, foreground="w")],
            color="#67000d",
            fontsize=8,
        )

        ax = available_regions.plot(
            projection=ccrs.PlateCarree(), add_ocean=True, text_kws=text_kws
        )

        plt.rcParams["figure.dpi"] = 200  # Set the DPI to 300 (adjust as needed)

        mask_plot = new_masks

        indices = export.where(export.tagging_region > 0, drop=True)

        ax.contourf(
            longitude, latitude, mask_plot, levels=[0.1, 1], zorder=5, alpha=0.8
        )

        if "prudence" in region_source:
            ax.set_extent([-20, 40, 30, 70])

        # Grid
        gl = ax.gridlines(draw_labels=True, alpha=0.5, linestyle=":", color="k")
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = {"size": fontsizes}
        gl.ylabel_style = {"size": fontsizes}

        return ax


def mask_with_shapefile(shapefile, regions, reference_file, output_dir, return_plot):
    """Mask tagging region using a shapefile.

    This function builds a tagging region based on an existing shapefile and stores the output file.

    Args:
        shapefile: name of the shapefile. Include: .dbf, .prj, .sbn, .sbx, .shp, .shp.xml and .shx files in the input directory.
        regions: if build up of multiple regions, indicate as list (int), otherwise FALSE.
        reference_file: a reference file of the preprocessing which is used to
            extract the dimensions of the data used.
        output_dir: path where to store the created tagging region.
        return_plot: if true, return a plot that shows the tagging region on a map.

    Returns:
        axes on which the region is plotted if 'return_plot' is True.
    """

    # Retrieve file dimensions
    file = xr.open_dataset(reference_file)

    longitude = np.array(file.longitude)
    latitude = np.array(file.latitude)

    # Read shapefile
    ds = gpd.read_file(shapefile)  # load shape file

    # Generate mask
    mask = regionmask.mask_3D_geopandas(ds, longitude, latitude)

    # Specify region ID - what to store
    if regions != False:
        region_ids = np.array(regions)
        region_ids = region_ids.tolist()

        if type(region_ids) == int:
            raise KeyError(
                "regions not specified correctly, please check whether your input is provided as a list [...,...], even if it is a single tagging region!"
            )

        new_masks = mask.isel(region=region_ids).values

        # Make WAM2layers fitting
        new_masks = np.where(
            new_masks == False, 0, 1
        )  # Replace True and False by 1 and 0's
        new_masks = np.nanmean(new_masks, axis=0)  # Multiple masks into 1D array
    else:
        new_masks = mask.isel(region=[0]).values

        # Make WAM2layers fitting
        new_masks = np.where(
            new_masks == False, 0, 1
        )  # Replace True and False by 1 and 0's
        new_masks = np.nanmean(new_masks, axis=0)  # Multiple masks into 1D array

    # Export as tagging region for WAM2layers

    # create xarray dataset
    data = new_masks

    export = xr.Dataset(
        {"tagging_region": (["latitude", "longitude"], data.astype(float))}
    )

    # set coordinates
    export["longitude"] = ("longitude", longitude)
    export["latitude"] = ("latitude", latitude)

    # save to NetCDF file

    export.to_netcdf(output_dir)

    logger.info(f"Stored tagging region in {output_dir}.")

    if return_plot:
        # Visualise all regions available

        fontsizes = 10
        pads = 20

        # Zoomed-in
        fig = plt.figure(figsize=(16, 10), dpi=200)
        ax = fig.add_subplot(111, projection=ccrs.PlateCarree())

        # Make figure
        ax.add_feature(cfeature.LAND, zorder=1, edgecolor="k", facecolor="w")

        ax.add_feature(cfeature.COASTLINE, zorder=2, linewidth=0.8)
        ax.add_feature(cfeature.BORDERS, zorder=2, linewidth=0.2, alpha=0.5)
        ax.add_feature(cfeature.RIVERS, zorder=2, linewidth=3, alpha=0.5)
        ax.add_feature(cfeature.STATES, zorder=2, facecolor="w")
        ax.add_feature(
            cfeature.LAKES,
            zorder=2,
            linewidth=0.8,
            edgecolor="k",
            alpha=0.5,
            facecolor="w",
        )
        mask_plot = new_masks

        indices = export.where(export.tagging_region > 0, drop=True)
        lon_indices = indices.longitude
        lon_indices_c = []
        for i in lon_indices:
            if i <= 180:
                item = i * -1
                lon_indices_c.append(item)
            else:
                lon_indices_c.append(i)

        ax.contourf(
            longitude, latitude, mask_plot, levels=[0.1, 1], zorder=5, alpha=0.8
        )
        ax.set_extent(
            [
                min(lon_indices_c) - 1,
                max(lon_indices_c) + 1,
                indices.latitude.min() - 1,
                indices.latitude.max() + 1,
            ]
        )
        ax.set_title("Zoomed-in", fontsize=fontsizes)

        # Grid
        gl = ax.gridlines(draw_labels=True, alpha=0.5, linestyle=":", color="k")
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = {"size": fontsizes * 2}
        gl.ylabel_style = {"size": fontsizes * 2}

        # Zoomed-out
        fig = plt.figure(figsize=(16, 10), dpi=200)
        ax2 = fig.add_subplot(111, projection=ccrs.PlateCarree())

        # Make figure
        ax2.add_feature(cfeature.LAND, zorder=1, edgecolor="k", facecolor="w")

        ax2.add_feature(cfeature.COASTLINE, zorder=2, linewidth=0.8)
        ax2.add_feature(cfeature.BORDERS, zorder=2, linewidth=0.2, alpha=0.5)
        ax2.add_feature(cfeature.RIVERS, zorder=2, linewidth=3, alpha=0.5)
        ax2.add_feature(cfeature.STATES, zorder=2, facecolor="w")
        ax2.add_feature(
            cfeature.LAKES,
            zorder=2,
            linewidth=0.8,
            edgecolor="k",
            alpha=0.5,
            facecolor="w",
        )
        mask_plot = new_masks

        indices = export.where(export.tagging_region > 0, drop=True)
        lon_indices = indices.longitude
        lon_indices_c = []
        for i in lon_indices:
            if i <= 180:
                item = i * -1
                lon_indices_c.append(item)
            else:
                lon_indices_c.append(i)

        ax2.contourf(
            longitude, latitude, mask_plot, levels=[0.1, 1], zorder=5, alpha=0.8
        )
        ax2.set_extent(
            [
                min(lon_indices_c) - 20,
                max(lon_indices_c) + 20,
                indices.latitude.min() - 20,
                indices.latitude.max() + 20,
            ]
        )
        ax2.set_title("Zoomed-out", fontsize=fontsizes)

        # Grid
        gl = ax2.gridlines(draw_labels=True, alpha=0.5, linestyle=":", color="k")
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = {"size": fontsizes * 2}
        gl.ylabel_style = {"size": fontsizes * 2}

        return ax, ax2


def mask_around_point(
    centerpoint,
    radius,
    reference_file,
    output_dir,
    return_plot,
):
    """Mask tagging region using a center point with a given radius.

    This function builds a square tagging region based on given coordinates and
    radius and stores it the output file.

    Args:
        centerpoint: Coordinates of the central point (latitude, longitude).
        radius (int): distance from center to edges of square box in degrees.
        reference_file: a reference file of the preprocessing which is used to
            extract the dimensions of the data used.
        output_dir: path where to store the created tagging region.
        return_plot: if true, return a plot that shows the tagging region on a map.

    Returns:
        axes on which the region is plotted if `return_plot' is True
    """

    # Retrieve file dimensions
    file = xr.open_dataset(reference_file)

    longitude_r = np.array(file.longitude)
    latitude_r = np.array(file.latitude)

    # Create a grid of latitudes and longitudes
    longitude, latitude = np.meshgrid(longitude_r, latitude_r)

    # retrieve latitude and longitude from ceterpoint
    center_lat, center_lon = centerpoint

    # Set the square region of value 1 in the mask
    mask = (np.abs(latitude - center_lat) <= radius) & (
        np.abs(longitude - center_lon) <= radius
    )

    new_masks = np.where(mask == False, 0, 1)  # Replace True and False by 1 and 0's

    # Export as tagging region for WAM2layers

    # create xarray dataset
    data = new_masks

    export = xr.Dataset(
        {"tagging_region": (["latitude", "longitude"], data.astype(float))}
    )

    # set coordinates
    export["longitude"] = ("longitude", longitude_r)
    export["latitude"] = ("latitude", latitude_r)

    # save to NetCDF file
    export.to_netcdf(output_dir)

    logger.info(f"Stored tagging region in {output_dir}.")

    if return_plot:
        # Visualise all regions available

        fontsizes = 10
        pads = 20

        fig = plt.figure(figsize=(16, 10), dpi=200)
        ax = fig.add_subplot(111, projection=ccrs.PlateCarree())

        # Make figure
        ax.add_feature(cfeature.LAND, zorder=1, edgecolor="k", facecolor="w")

        ax.add_feature(cfeature.COASTLINE, zorder=2, linewidth=0.8)
        ax.add_feature(cfeature.BORDERS, zorder=2, linewidth=0.2, alpha=0.5)
        ax.add_feature(cfeature.RIVERS, zorder=2, linewidth=3, alpha=0.5)
        ax.add_feature(cfeature.STATES, zorder=2, facecolor="w")
        ax.add_feature(
            cfeature.LAKES,
            zorder=2,
            linewidth=0.8,
            edgecolor="k",
            alpha=0.5,
            facecolor="w",
        )
        mask_plot = new_masks

        ax.contourf(longitude_r, latitude_r, mask_plot, levels=[0.1, 1])

        indices = export.where(export.tagging_region > 0, drop=True)
        ax.contourf(
            longitude, latitude, mask_plot, levels=[0.1, 1], zorder=5, alpha=0.8
        )
        ax.set_extent(
            [
                indices.longitude.min() - 10,
                indices.longitude.max() + 10,
                indices.latitude.min() - 10,
                indices.latitude.max() + 10,
            ]
        )

        ax.plot(
            center_lon,
            center_lat,
            marker="o",
            markerfacecolor="r",
            markeredgecolor="k",
            linestyle="None",
            markersize=radius * 0.5,
            zorder=6,
        )

        # Grid
        gl = ax.gridlines(draw_labels=True, alpha=0.5, linestyle=":", color="k")
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = {"size": fontsizes * 2}
        gl.ylabel_style = {"size": fontsizes * 2}
        return ax


def mask_around_track(
    centerpoints,
    times,
    radius,
    reference_file,
    output_file,
    return_plot,
):
    """Mask tagging region using a track of center points with a given radius (currently in degrees).

    This function builds a square source region based on given coordinates and
    radius and stores it the output file.

    Args:
        centerpoints: List of coordinates of the central points, changing over time [(lat1, lon1), ....,lat_final, lon_final)] .
        radius (int): distance from center to edges of square box in degrees.
        reference_file: a reference file of the preprocessing which is used to
            extract the dimensions of the data used.
        output_file: path where to store the created source region.
        return_plot: if true, return a plot that shows the source region on a map.

    Returns:
        axes on which the region is plotted if `return_plot' is True
    """

    # Retrieve file dimensions
    file = xr.open_dataset(reference_file)

    longitude_r = np.array(file.longitude)
    latitude_r = np.array(file.latitude)

    # Create a grid of latitudes and longitudes
    longitude, latitude = np.meshgrid(longitude_r, latitude_r)

    # Initiate result array
    new_masks = np.zeros((len(centerpoints), len(latitude_r), len(longitude_r)))

    # loop over times
    for t in range(len(centerpoints)):
        centerpoint = centerpoints[t]
        # retrieve latitude and longitude from ceterpoint
        center_lat, center_lon = centerpoint

        # Set the square region of value 1 in the mask
        mask = (np.abs(latitude - center_lat) <= radius) & (
            np.abs(longitude - center_lon) <= radius
        )

        new_masks[t, ::] = np.where(
            mask == False, 0, 1
        )  # Replace True and False by 1 and 0's

    # Export as source region for WAM2layers
    # create xarray dataset
    data = new_masks

    export = xr.Dataset(
        {"source_region": (["time", "latitude", "longitude"], data.astype(float))}
    )

    # set coordinates
    export["longitude"] = ("longitude", longitude_r)
    export["latitude"] = ("latitude", latitude_r)
    export["time"] = ("time", times)

    # save to NetCDF file
    export.to_netcdf(output_file)

    logger.info(f"Stored source region in {output_file}.")

    if return_plot:
        # Visualise all regions available

        fontsizes = 10
        pads = 20

        fig = plt.figure(figsize=(16, 10), dpi=200)
        ax = fig.add_subplot(111, projection=ccrs.PlateCarree())

        # Make figure
        ax.add_feature(cfeature.LAND, zorder=1, edgecolor="k", facecolor="w")

        ax.add_feature(cfeature.COASTLINE, zorder=2, linewidth=0.8)
        ax.add_feature(cfeature.BORDERS, zorder=2, linewidth=0.2, alpha=0.5)
        ax.add_feature(cfeature.RIVERS, zorder=2, linewidth=3, alpha=0.5)
        ax.add_feature(cfeature.STATES, zorder=2, facecolor="w")
        ax.add_feature(
            cfeature.LAKES,
            zorder=2,
            linewidth=0.8,
            edgecolor="k",
            alpha=0.5,
            facecolor="w",
        )
        mask_plot = new_masks

        ax.contourf(longitude_r, latitude_r, mask_plot, levels=[0.1, 1])

        indices = export.where(export.source_region == 1, drop=True)
        ax.contourf(
            longitude, latitude, mask_plot, levels=[0.1, 1], zorder=5, alpha=0.8
        )
        ax.set_extent(
            [
                indices.longitude.min() - 10,
                indices.longitude.max() + 10,
                indices.latitude.min() - 10,
                indices.latitude.max() + 10,
            ]
        )

        ax.plot(
            center_lon,
            center_lat,
            marker="o",
            markerfacecolor="r",
            markeredgecolor="k",
            linestyle="None",
            markersize=radius * 0.5,
            zorder=6,
        )

        # Grid
        gl = ax.gridlines(draw_labels=True, alpha=0.5, linestyle=":", color="k")
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = {"size": fontsizes * 2}
        gl.ylabel_style = {"size": fontsizes * 2}
        return ax
