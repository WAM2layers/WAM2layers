"""
WAM2-layers Source Region Desinger

@author: Vincent de Feiter (vincent.defeiter@wur.nl) | GitHub: vincentdefeiter

Version 2.0 | 26-03-2023

Welcome to the WAM2-Layers Source Region Designer. With this designer, a source region for the WAM2-Layers model will be made. A source region,
is defined as a region with zeros and ones. The regions with the value 1 are recognized by the model to track the precipitation events from this region solely.
The regions where from where no precipitation should be tracked, are denoted as 0.

In this designer, a total of 3 methods (functions) can be used. See under FUNCTIONS if you would like to change the contents of the functions.
Otherwise, you can directly go to the RUN codelines.

Do not forget to provide your OUTPUT directory just below in GENERAL settings!

The methods are:

    > Method 1: Using the available regions from the 'region mask' python package.

    Here, available regions from the 'Region Mask' Python Package (https://regionmask.readthedocs.io/en/stable/defined_scientific.html)
    are used. You are able to specify 1 or multiple source regions.

    > Method 2: Using a shapefile combined with 'geopandas' and the 'region mask' python package.

    Here, an available source region shapefile (downloaded from an online source) can be implemented.
    You are avaialble to specify 1 specific region, or connect multiple (or multiple layers) of the
    shapefile as source region.

    Shapefiles can be retrieved from various sources, a good source is the HydroSHEDS data of the World Wide
    Fund for Nature (WWF). The HydroSHEDS data supplies shapefiles of numerous river basins:
    Can be obtained via https://www.hydrosheds.org/products/hydrobasins.

    Lehner, B., Grill G. (2013): Global river hydrography and network routing: baseline data and
    new approaches to study the world’s large river systems. Hydrological Processes, 27(15):
    2171–2186. Data is available at www.hydrosheds.org.

    NOTE: The data from HydroSHEDS is provided in Pfafstetter coded level. The shapefile becomes more complex,
    e.g., more smaller regions are drawn, with each increment in level (higher level). See documentation for
    more information.

    > Method 3: Center point and a squared perimeter around it.

    Here, a set center point (see settings) is set, after which a squared region (depending on the number of degrees)
    is drawn around this center point. Take care that no emphasis is set on borders between sea/land.
"""

import cartopy.crs as ccrs  # used for plotting on a map
import geopandas as gpd  # reading shapefiles
import matplotlib.patheffects as pe  # additional plotting functionalities
import matplotlib.pyplot as plt  # used for plotting purposes
import numpy as np  # used for array calculations
import regionmask  # used for region mapping
import xarray as xr  # used to read .nc files
from cartopy import feature as cfeature


def mask_with_regionmask(
    region_source, region_names, preprocessfile, output_directory, plotting
):
    """
    Source region by using a named region from regionmask.

    Here, available regions from the 'Region Mask' Python Package
    (https://regionmask.readthedocs.io/en/stable/defined_scientific.html)
    are used. You are able to specify 1 or multiple source regions, indicated with a list [...].

    The properties in this function indicate the following:

    region_source:
    Please indicate here a region source from the Region Mask Python Package (as a STRING!), which can be:
            > Giorig Regions: The Giorgi reference regions, rectangular regions proposed in Giorgi and Francisco, 2000 were used in the third and fourth assessment reports of the Intergovernmental Panel on Climate Change (IPCC).
            > SREX Regions: The Giorgi regions were modified using more flexible polygons in the IPCC Special Report on Managing the Risks of Extreme Events and Disasters to Advance Climate Adaptation (SREX; Seneviratne et al., 2012). A slightly modified and extended version of the SREX regions (not included in regionmask) was used for the fifth IPCC assessment report.
            > AR6 Regions: The sixth IPCC assessment report (AR6) again updated the SREX regions in Iturbide et al., (2020), which defines 58 regions. The regions cover the land and ocean (ar6.all). In addition the regions are also divided into land (ar6.land) and ocean (ar6.ocean) categories. The numbering is kept consistent between the categories. Note that some regions are in the land and in the ocean categories (e.g. the Mediterranean).
            > PRUDENCE Regions: The PRUDENCE regions were defined in the PRUDENCE project as European sub-areas for regional climate model output and are often used in European climate studies. They contain 8 regions, Alps (AL), British Isles (BI), Eastern Europe (EA), France (FR), Iberian Peninsula (IP), Mediterranean (MD), Mid-Europe (ME), and Scandinavia (SC).

            For more information on which specific regions, check the documentation of the Region Mask Python Package.

    ###########################################################################

    region_names:
    Please provide here a list '[...]' with the value(s) of the regions of interest. For example, [50,31,32]. See the package information for more details.

    ###########################################################################

    preprocessfile:

    Please provide here a reference file of the preprocessing which is used to extract the dimensions of the data used. This step is required to put the source region back into the data dimensions.

    ###########################################################################

    plotting:
    This option can be turned 'True' or 'False' depending on whether you want to plot your source region to see the end result. This can be used as a visual inspection to verify that the correct soruce region is generated.

    """

    # Retrieve file dimensions
    file = xr.open_dataset(preprocessfile)

    longitude = np.array(file.longitude)
    latitude = np.array(file.latitude)

    # Set available regions
    try:
        available_regions = getattr(regionmask.defined_regions, region_source)
    except:
        raise KeyError(
            'region_source not specified correctly, please check the "regionmask.defined_regions" options in the Region Mask documentation for more details'
        )

    # create (a) mask(s) of the required region(s):
    mask = getattr(available_regions, "mask_3D")(longitude, latitude)

    # Specify region ID - what to store
    region_ids = np.array(region_names) - 1
    region_ids = region_ids.tolist()

    if type(region_ids) == int:
        raise KeyError(
            "region_names not specified correctly, please check whether your input is provided as a list [...,...], even if it is a single source region!"
        )

    new_masks = mask.isel(region=region_ids).values

    # Make WAM2Layers fitting
    new_masks = np.where(
        new_masks == False, 0, 1
    )  # Replace True and False by 1 and 0's
    new_masks = np.nanmean(new_masks, axis=0)  # Multiple masks into 1D array

    # Plot as a check
    if plotting:
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

        mask_plot = new_masks

        ax.contourf(longitude, latitude, mask_plot, levels=[0.1, 1])
        ax.set_ylim(latitude[0], latitude[-1])
        ax.set_xlim(longitude[0], longitude[-1])
        ax.invert_yaxis()

        # Grid
        gl = ax.gridlines(draw_labels=True, alpha=0.5, linestyle=":", color="k")
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = {"size": fontsizes}
        gl.ylabel_style = {"size": fontsizes}

        plt.savefig(
            output_directory + "mask_with_regionmask_selected_region(s).png",
            dpi=200,
            bbox_inches="tight",
        )

    # Export as source region for WAM2Layers

    # create xarray dataset
    data = new_masks

    export = xr.Dataset(
        {"source_region": (["latitude", "longitude"], data.astype(float))}
    )

    # set coordinates
    export["longitude"] = ("longitude", longitude)
    export["latitude"] = ("latitude", latitude)

    # save to NetCDF file
    export.to_netcdf(output_directory + "source_region.nc")
    return ax


def mask_with_shapefile(
    shapefile, region_names, preprocessfile, output_directory, plotting
):
    """
    Mask source region using a shapefile.

    This function builds a source region based on an existing shapefile. It will make use of the 'Region Mask' Python Package
    (https://regionmask.readthedocs.io/en/stable/defined_scientific.html).

    The properties in this function indicate the following:

    shapefile:
    Please indicate here the source of your shapefile. Make sure that you also have the .dbf, .prj, .sbn
    .sbx, .shp, .shp.xml and .shx files stored in the same directory as your shapefile.

    ###########################################################################

    region_names:
    Some shapefiles are build up of different regions, e.g., provinces or basins. Please provide here a list '[...]' with the value(s) of the regions of interest. For example, [50,31,32]. See the shapefile documentation for more information.
    If there are no different regions in the shapefile, please specifiy this property as: FALSE.

    ###########################################################################

    preprocessfile:

    Please provide here a reference file of the preprocessing which is used to extract the dimensions of the data used. This step is required to put the source region back into the data dimensions.

    ###########################################################################

    plotting:
    This option can be turned 'True' or 'False' depending on whether you want to plot your source region to see the end result. This can be used as a visual inspection to verify that the correct soruce region is generated.

    """
    # Retrieve file dimensions
    file = xr.open_dataset(preprocessfile)

    longitude = np.array(file.longitude)
    latitude = np.array(file.latitude)

    # Read shapefile
    ds = gpd.read_file(shapefile)  # load shape file

    # Generate mask
    mask = regionmask.mask_3D_geopandas(ds, longitude, latitude)

    # Specify region ID - what to store
    if region_names != False:
        region_ids = np.array(region_names)
        region_ids = region_ids.tolist()

        if type(region_ids) == int:
            raise KeyError(
                "region_names not specified correctly, please check whether your input is provided as a list [...,...], even if it is a single source region!"
            )

        new_masks = mask.isel(region=region_ids).values

        # Make WAM2Layers fitting
        new_masks = np.where(
            new_masks == False, 0, 1
        )  # Replace True and False by 1 and 0's
        new_masks = np.nanmean(new_masks, axis=0)  # Multiple masks into 1D array
    else:
        new_masks = mask.isel(region=[0]).values

        # Make WAM2Layers fitting
        new_masks = np.where(
            new_masks == False, 0, 1
        )  # Replace True and False by 1 and 0's
        new_masks = np.nanmean(new_masks, axis=0)  # Multiple masks into 1D array

    # Plot as a check
    if plotting == True:
        # Visualise all regions available

        fontsizes = 10
        pads = 20

        fig = plt.figure(figsize=(16, 10))
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

        ax.contourf(longitude, latitude, mask_plot, levels=[0.1, 1])
        ax.set_ylim(latitude[0], latitude[-1])
        ax.set_xlim(longitude[0], longitude[-1])
        ax.invert_yaxis()

        # Grid
        gl = ax.gridlines(draw_labels=True, alpha=0.5, linestyle=":", color="k")
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = {"size": fontsizes}
        gl.ylabel_style = {"size": fontsizes}

        plt.savefig(
            output_directory + "mask_with_shapefile_selected_region(s).png",
            dpi=200,
            bbox_inches="tight",
        )

    # Export as source region for WAM2Layers

    # create xarray dataset
    data = new_masks

    export = xr.Dataset(
        {"source_region": (["latitude", "longitude"], data.astype(float))}
    )

    # set coordinates
    export["longitude"] = ("longitude", longitude)
    export["latitude"] = ("latitude", latitude)

    # save to NetCDF file
    export.to_netcdf(output_directory + "source_region.nc")

    return ax


def mask_around_point(
    centerpoint,
    radius,
    reference_file,
    output_file="source_region.nc",
    return_plot=False,
):
    """Mask source region using a center point with a given radius.

    This function builds a square source region based on given coordinates and
    radius and stores it the output file.

    Args:
        centerpoint: Coordinates of the central point (latitude, longitude).
        radius (int): distance from center to edges of square box in degrees.
        reference_file: a reference file of the preprocessing which is used to
            extract the dimensions of the data used.
        output_file: path where to store the regionmask file.
        return plot: if true, return a plot that shows the source region on a map.

    Returns:
        axes on which the region is plotted if `plotted` is True
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

    # Export as source region for WAM2Layers

    # create xarray dataset
    data = new_masks

    export = xr.Dataset(
        {"source_region": (["latitude", "longitude"], data.astype(float))}
    )

    # set coordinates
    export["longitude"] = ("longitude", longitude_r)
    export["latitude"] = ("latitude", latitude_r)

    # save to NetCDF file
    export.to_netcdf(output_file)

    print(f"Stored source region in {output_file}.")

    if return_plot:
        # Visualise all regions available

        fontsizes = 10
        pads = 20

        fig = plt.figure(figsize=(16, 10))
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
        ax.set_ylim(latitude_r[0], latitude_r[-1])
        ax.set_xlim(longitude_r[0], longitude_r[-1])
        ax.invert_yaxis()

        ax.plot(
            center_lon,
            center_lat,
            marker="o",
            markerfacecolor="r",
            markeredgecolor="k",
            linestyle="None",
            markersize=radius * 0.5,
        )

        # Grid
        gl = ax.gridlines(draw_labels=True, alpha=0.5, linestyle=":", color="k")
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = {"size": fontsizes}
        gl.ylabel_style = {"size": fontsizes}
        return ax
