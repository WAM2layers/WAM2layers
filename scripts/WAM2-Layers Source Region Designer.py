"""
WAM2-layers Source Region Desinger

@author: Vincent de Feiter (vincent.defeiter@wur.nl)

version 1.0 | 19-02-2023

"""
#%%
#---------------------------------------------------------------------------
#               I  M  P  O  R  T   M  O D U L E S 
#---------------------------------------------------------------------------
#%%
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import regionmask
import matplotlib.patheffects as pe
import geopandas as gpd


#%%
#---------------------------------------------------------------------------
#              S O U R C E  R E G I O N  D E S I G N E R S 
#---------------------------------------------------------------------------
#%%

'''
Welcome to the WAM2-Layers Source Region Designer. All methods are dependent on the 
resolution of your ERA5 (or other) data, specified in the general settings. 

There are 3 methods specified below:
    
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
    

To set the specific method, see the settings below. Please note that you are able to specify more of the settings within the original code. 

'''

'''

S E T T I N G S 

'''

#=======================================
            #GENERAL SETTINGS
#=======================================

#Output directory - Set here where the output is stored
output_dir = '...'

#Specify latitude and longitude - same as input data to WAM2Layers
lat_max = 80
lat_min = -80

lon_max = 180
lon_min = -180

step = 0.25 #0.25 degree intervals (ERA5 data)

#Ranges
latitude = np.arange(lat_min,lat_max+step,step) #check ranges based on your input
longitude = np.arange(lon_min,lon_max,step)

#=======================================
            #METHOD SELECTOR
#=======================================
selector = 1    #1 = Method 1, 2 = Method 2, 3 = Method 3. 

#======================================================================================================================

#=======================================
        #METHOD SPECIFIC SETTINGS
#=======================================

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#           M E T H O D  1
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#Create an overview of available regions from: https://regionmask.readthedocs.io/en/stable/defined_scientific.html
available_regions = regionmask.defined_regions.giorgi #<-- specify available regions source (giorgi or other)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#           M E T H O D  2
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
shapefile = '...'    #specify a shapefile

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#           M E T H O D  3
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Define the center point and distance in degrees

#Center point
center_lat = -2         
center_lon = -58.75
    
distance = 20               #Distance of the square around the center point, in degrees


'''

PLEASE NOTE: region ID's (what specific region or regions) to export as a source region, need to be specified within the code itself. 
A visual representation of the region is generated to show intermediate steps / the final product. 

'''


#%%

'''

METHOD 1 - Using the available regions from the 'region mask' python package
 
'''

if selector == 1:
        
    #Visualisation of the regions and their ID's --> select region
    text_kws = dict(
        bbox=dict(color="none"),
        path_effects=[pe.withStroke(linewidth=2, foreground="w")],
        color="#67000d",
        fontsize=8,
    )
    
    ax = available_regions.plot(
        projection=ccrs.Robinson(), add_ocean=True, text_kws=text_kws
    )
    
    ax.set_global()   
    plt.show()
    
    
    #create (a) mask(s) of the required region(s):
    mask = regionmask.defined_regions.giorgi.mask_3D(longitude, latitude)   #<-- specify available regions source
    
    #Specify region ID - what to store
    region_id = 2                                                           #single mask id
    region_ids = [1]                                                        #multiple regions in one mask [id-1,id-1,...]
    ''''!! NOTE REGION_IDS is INDEX of map - 1 !! '''
    
    new_mask = mask.sel(region=region_id).values 
    new_masks = mask.isel(region=region_ids).values
    
    #Make WAM2Layers fitting
    new_mask = np.where(new_mask == False, 0, 1)                            #Replace True and False by 1 and 0's
    new_masks = np.where(new_masks == False, 0, 1)                          #Replace True and False by 1 and 0's
    new_masks = np.nanmean(new_masks,axis=0)                                #Multiple masks into 1D array
    
    #Plot as a check
    plt.contourf(longitude,latitude,new_masks)                              #check (switch between mask or masks when needed)
    plt.show()
    
    #Export as source region for WAM2Layers
    
    # create xarray dataset
    data = new_mask   #switch between mask or masks
    
    export = xr.Dataset({'source_region': (['latitude', 'longitude'], data.astype(float))})
        
    # set coordinates
    export['longitude'] = ('longitude', longitude)
    export['latitude'] = ('latitude', latitude)
    
    # save to NetCDF file
    export.to_netcdf(output_dir+'source_region.nc')
    

#%%

'''

METHOD 2 - Using a shapefile combined with 'geopandas' and the 'region mask' python package
 
'''

if selector == 2:
    ds = gpd.read_file(shapefile)                    #load shape file
    
    #visualisation
    ds.boundary.plot(color=None, edgecolor="k", linewidth=0.8)      #check
    plt.show()
    
    #Generate mask
    mask = regionmask.mask_3D_geopandas(ds, longitude, latitude)
    
    #Specify region ID - what to store
    region_id = 1                                                         #single mask id
    region_ids = [1]                                                      #multiple regions in one mask [id-1,id-1,...]
    ''''!! NOTE REGION_IDS is INDEX of map - 1 !! '''
    
    new_mask = mask.sel(region=region_id).values 
    new_masks = mask.isel(region=region_ids).values
    
    #Make WAM2Layers fitting
    new_mask = np.where(new_mask == False, 0, 1)                            #Replace True and False by 1 and 0's
    new_masks = np.where(new_masks == False, 0, 1)                          #Replace True and False by 1 and 0's
    new_masks = np.nanmean(new_masks,axis=0)                                #Multiple masks into 1D array
        
    #Plot as a check
    plt.contourf(longitude,latitude,new_mask)                               #check (switch between mask or masks when needed)
    plt.show()
    
    #Export as source region for WAM2Layers
    
    # create xarray dataset
    data = new_mask   #switch between mask or masks
    
    export = xr.Dataset({'source_region': (['latitude', 'longitude'], data.astype(float))})
        
    # set coordinates
    export['longitude'] = ('longitude', longitude)
    export['latitude'] = ('latitude', latitude)
    
    # save to NetCDF file
    export.to_netcdf(output_dir+'source_region.nc')




#%%

'''

METHOD 3 - Center point and a squared perimeter around it
 
'''
        
if selector == 3:
        
    latitude_r = latitude
    longitude_r = longitude
    
    # Create a grid of latitudes and longitudes
    longitude, latitude = np.meshgrid(longitude,latitude)
    
    # Set the square region of value 1 in the mask
    mask = (np.abs(latitude - center_lat) <= distance) & (np.abs(longitude - center_lon) <= distance)
    
    new_mask = np.where(mask == False, 0, 1)                      #Replace True and False by 1 and 0's
    
    #plot
    plt.contourf(longitude,latitude,new_mask)                     #check - Note dimensions of the axis may not match!
    plt.plot(center_lon,center_lat,'o')
    plt.show()
    
    #Export as source region for WAM2Layers
    
    # create xarray dataset
    data = new_mask
    
    export = xr.Dataset({'source_region': (['latitude', 'longitude'], data.astype(float))})
        
    # set coordinates
    export['longitude'] = ('longitude', longitude_r)
    export['latitude'] = ('latitude', latitude_r)
    
    # save to NetCDF file
    export.to_netcdf(output_dir+'source_region.nc')
