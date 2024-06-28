# coding: cp1252
''' 

Visualisation Functions WAM-2Layers

Version 2.0 | 14-07-2023

Updated on 17-06-2024

@author: Vincent de Feiter (vincent.defeiter@wur.nl) | GitHub: vincentdefeiter

Welcome to this visualisation script. It provides an overview of functions that can
be utilised to plot different input or output data and model results. The functions 
serve as a 'quick look' or general idea/basis to continue your own development. 

The following functions are prescribed: 

> Logfile. This function provides a logfile (.txt) with the evaporation recycling ratio (evaporated water that returns as precipitation), pecipitation recycling ratio (region�s dependence on evaporation from within the region for precipitation), sum of the precipitation, evaporation and their respective absolute differnce. This distance should be kept to a minimum and can be used as a check whether the model is operated correctly. 

> Precipitation. This function develops a figure showing the precipitation events tracked by WAM2Layers. It visualises the source region accompanied by the precipitation. 

> Evaporation. This function develops a figure showing the tracked evaporation (showing source and sinks) from the source region, with the accompanying moisture fluxes (arrows) and mean surface pressure (contours in hPa).

> Moisture_fluxes. This function develops figures showing the mean moisture fluxes (arrows) with their corresponding magnitude (contours). It develops a figure for a specific time range and for both the lower and upper layer of the model (seperated by the p_divide). 



'''
from pathlib import Path
import click
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import AutoMinorLocator
import xarray as xr
import pandas as pd
from cartopy import crs
from cartopy import feature as cfeature
from cmocean import cm
from wam2layers.config import Config
from wam2layers.utils.grid import get_grid_info
from wam2layers.tracking.io import input_path, load_tagging_region, output_path
from scipy.interpolate import griddata
import netCDF4 as nc
import numpy as np
import xarray as xr
import glob
import os
import datetime
import math
import cartopy.mpl.ticker as cticker
import matplotlib.patches as mpatches
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta


def logfile(config_file,filename,outdir):
    
    '''
    Logfile (.txt) with the evaporation recycling ratio, pecipitation recycling ratio, sum of the precipitation, evaporation and their respective absolute differnce. 
    
    Args:
      config_file (str): path to the .yaml configuration file to run the model. 
      filename (str): string to specify the outputfile of the configfile on XXX --> 'Logfile_XXXX.txt'
      outdir (str): path to an output directory to store the output config file. 
      
    Returns:
      Logfile (.txt) with the requested information.
    
    '''
        
    #Load general data
    config = Config.from_yaml(config_file)
    region = load_tagging_region(config)
    
    #Load tracked moisture and precipitation
    a_gridcell, lx, ly = get_grid_info(region)
    output_files = f"{config.output_folder}/*.nc"
    ds_e = xr.open_mfdataset(output_files, combine='nested', concat_dim='time')
    e_track = ds_e.e_track.sum('time').compute()
    p_tagged = ds_e.tagged_precip.sum('time').compute()

    #calculate precipitation recycling ratio + mean precipitation over area
    p_r=(e_track*a_gridcell[:, None]).where(region, drop=True).sum() / (p_tagged*a_gridcell[:, None]).sum()  
 
    #Calculate (area-weighted) mean precip and moisture sources 
    precip_mean  = (p_tagged*a_gridcell[:, None]).sum()/a_gridcell[:, None].sum() #Should only be target area
    e_track_mean  = (e_track*a_gridcell[:, None]).sum()/a_gridcell[:, None].sum()
    
    diff_areacorrected = (p_tagged*a_gridcell[:, None]).sum() - (e_track*a_gridcell[:, None]).sum()
    
    #Write to file
    logfile = open(outdir+"/Logfile_"+filename+".txt",'w')
    logfile.write('mean precipiation: '+str(float(precip_mean))+'\n')
    logfile.write('mean moisture sources: '+str(float(e_track_mean))+'\n')
    logfile.write('Difference: '+str(float(diff_areacorrected))+'\n')
    logfile.write('p_r: '+str(float(p_r*100))+' % \n')
    
    logfile.close()
    

def precipitation(config_file,filename,outdir,step, fontsizes, pads):
    """
    
    A figure showing the precipitation events tracked by WAM-2Layers
    
    Args:
      config_file (str): path to the .yaml configuration file to run the model. 
      filename (str): string to specify the output file.'
      outdir (str): path to an output directory to store the output config file. 
      step (int): stepsize of the ticks for the colourbar. 
      fontsizes (int): fontsizes. 
      pads (int): integer indicating the spacing between the labels and the axes. 
      
    Returns:
      figure (.png) showing the precipitation events tracked by WAM-2Layers.
        
    
    """
    
    #Load data
    config = Config.from_yaml(config_file)
    region = load_tagging_region(config)
    
    input_files = f"{config.preprocessed_data_folder}/*.nc"
    ds = xr.open_mfdataset(input_files, combine='nested', concat_dim='time')
    start = config.tagging_start_date
    end = config.tagging_end_date
        
    #Load precipitaiton
    subset = ds.precip.sel(time=slice(start, end))
    precip = (subset * region * 3600).sum('time').compute()
    precip = np.where(precip>1,precip,np.nan)
    
    #Get limits
    longitude = np.array(region.longitude.values)
    latitude = np.array(region.latitude.values)

    #Design colourmap
    cmap_name = 'my_list'
    colors = ['w','#00008B','#800080','#E50000','#FFD700']  
    bin = 50
    cmap_own = LinearSegmentedColormap.from_list(cmap_name, colors, N=bin)
    
    # Make figure
    fig = plt.figure(figsize=(16, 10))
    ax = fig.add_subplot(111, projection=crs.PlateCarree())

    # Add cartopy information
    ax.add_feature(cfeature.LAND, zorder=1, edgecolor='k',facecolor='w')
    ax.add_feature(cfeature.OCEAN, zorder=1, facecolor='k',alpha=0.1)
    ax.add_feature(cfeature.COASTLINE,zorder=2,linewidth=0.8)
    ax.add_feature(cfeature.BORDERS,zorder=2,linewidth=0.2,alpha=0.5)
    ax.add_feature(cfeature.RIVERS,zorder=2,linewidth=3,alpha=0.5)
    ax.add_feature(cfeature.STATES,zorder=2,facecolor='w')
    ax.add_feature(cfeature.LAKES,zorder=2,linewidth=0.8,edgecolor='k',alpha=0.5,facecolor='w')
    
    #Plot precipitation
    level_max = (np.round(np.nanmax(precip))//step)*step
    level_min = 0
    levels_cb = np.arange(0,level_max+1,1)
    levels_ticks = np.arange(0,level_max+step,step)
    
    cb = ax.contourf(ds.coords['longitude'],ds.coords['latitude'],precip,cmap=cmap_own, levels =levels_cb,zorder=3)
    cbar = plt.colorbar(cb, ticks=levels_ticks)
    cbar.ax.set_title('(mm)', fontsize=fontsizes,pad=pads*0.5)
    cbar.ax.tick_params(labelsize=fontsizes)    
    
    #Add border of the configfile
    region.plot.contour(ax=ax, levels = [0.1,1], linewidths=3,colors='k',zorder=4)
    
    #Plot limit  
    index = np.where(region.values==1)
    latitude_sel = latitude[index[0]]
    longitude_sel = longitude[index[1]]
      
    ax.set_ylim(np.nanmin(latitude_sel)-5,np.nanmax(latitude_sel)+5)
    
    converted_longitudes = (longitude_sel + 180) % 360 - 180
    ax.set_xlim(np.nanmin(converted_longitudes)-5,np.nanmax(converted_longitudes)+5)
       
    #Add minimap
    cmap_name = 'my_list'
    colors = ['tab:blue','tab:blue']  
    bin = 50
    cmap_own = LinearSegmentedColormap.from_list(cmap_name, colors, N=bin)
    
    #ax_mini = fig.add_axes([0.01, 0.01, 0.12, 0.2], projection=crs.PlateCarree(),facecolor='None') #x,y,w,h
    ax_mini = fig.add_axes([0.1, 0.1, 0.12, 0.2], projection=crs.PlateCarree(),facecolor='None') #x,y,w,h

    ax_mini.set_extent([np.nanmin(converted_longitudes)-50,np.nanmax(converted_longitudes)+50, np.nanmin(latitude)-50,np.nanmax(latitude)+50], crs=crs.PlateCarree())
    ax_mini.coastlines(alpha=0.5)
    ax_mini.add_feature(cfeature.COASTLINE,alpha=0.5)
    ax_mini.add_feature(cfeature.BORDERS,alpha=0.5)
    ax_mini.contourf(longitude,latitude,region.values, levels = [0.1,1])
    ax_mini.contour(longitude,latitude,region.values, levels = [0.1,1], linewidths=3,colors='k')
    
    #Grid
    gl = ax.gridlines(draw_labels=True,alpha=0.5, linestyle=':', color='k')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': fontsizes}
    gl.ylabel_style = {'size': fontsizes}
    
    # Export figure
    out_dir = Path(config.output_folder) / "figures"
    out_dir.mkdir(exist_ok=True, parents=True)
    fig.savefig(outdir+"/Precipitation_input_"+filename+".png", bbox_inches='tight', dpi=200)

def evaporation(config_file,filename,data_path,outdir,step, fontsizes, pads):
    """
    
    A figure showing the evaporation tracked by WAM-2Layers
    
    Args:
      config_file (str): path to the .yaml configuration file to run the model. 
      filename (str): string to specify the output file.
      data_path (str): path to the location where the ERA5 (or other) data is stored. 
      Make sure to have data per year and per month. 
      
      outdir (str): path to an output directory to store the output config file. 
      step (int): stepsize of the ticks for the colourbar. 
      fontsizes (int): fontsizes. 
      pads (int): integer indicating the spacing between the labels and the axes. 
      
    Returns:
      figure (.png) showing the evaporation tracked by WAM-2Layers.
        
    
    """

    #Load data
    config = Config.from_yaml(config_file)
    region = load_tagging_region(config)
    a_gridcell, lx, ly = get_grid_info(region)
    output_files = f"{config.output_folder}/*.nc"
    ds = xr.open_mfdataset(output_files, combine='nested', concat_dim='time')
    
    start = config.tracking_start_date
    end = config.tracking_end_date
    
    #Load evaporation
    e_track = ds.e_track.sum('time').compute()
    e_track = np.where(e_track>0.1,e_track,np.nan)
        
    #Background data - isobars, wind field etc. 
   
    #Start date
    start_date = str(start)
    start_date = datetime.strptime(start_date, '%Y-%m-%d %H:%M:%S')
    
    #End date 
    end_date = str(end)
    end_date = datetime.strptime(end_date, '%Y-%m-%d %H:%M:%S')
            
    date_time = []
    one_month_delta = relativedelta(months=1)
    current_date = start_date
    while current_date <= end_date:
        string = current_date.strftime('%Y%m%d')
        date_time.append(datetime.strptime(string,'%Y%m%d'))
        current_date += one_month_delta #timedelta(month=1)
    
    months=[]
    years=[]
    for item in date_time:
        months.append(str(item)[5:7])
        years.append(str(item)[:4])
    
    
    #Surface pressure
    print('now loading p')
  
    all_sp = []
    for year in years:
        for month in months:
            ERA5_data_path = data_path+year+'/'+month+'/'
            subs = '*sp.nc'
            ds = xr.open_mfdataset(ERA5_data_path+subs)
            ds = ds / 100
            all_sp.append(ds)

    ds_sp = xr.concat(all_sp, dim="time")
    mean = ds_sp.mean(dim='time')
    
    print('p loaded')
    
    #Moisture fluxes   
    path_input = f"{config.preprocessed_data_folder}/*.nc"
    ds = xr.open_mfdataset(path_input)
    subset = ds.sel(time=slice(start, end))
      
    mean_subset = subset.mean(dim='time')        
    mean_subset = mean_subset.assign(mean_lower=(np.sqrt((mean_subset['fy_lower'])**2 + (mean_subset['fx_lower'])**2)))    
    mean_subset = mean_subset.assign(mean_upper=(np.sqrt((mean_subset['fy_upper'])**2 + (mean_subset['fx_upper'])**2)))    
    mean_subset = mean_subset.assign(mean_x=(mean_subset['fx_upper'] + mean_subset['fx_lower']))    
    mean_subset = mean_subset.assign(mean_y=(mean_subset['fy_upper'] + mean_subset['fy_lower']))
    
    #Get limits
    longitude = np.array(region.longitude.values)
    latitude = np.array(region.latitude.values)
    
    #own colourmap
    cmap_name = 'my_list'
    colors = ['w','#00008B','#800080','#E50000','#FFD700']  
    bin = 200
    cmap_own = LinearSegmentedColormap.from_list(cmap_name, colors, N=bin)

    # Make figure
    fig = plt.figure(figsize=(16, 10))
    ax = fig.add_subplot(111, projection=crs.PlateCarree())
    
    #pressure
    min_p = (np.round(np.nanmean(mean['sp']))//5)*5
    max_p = (np.round(np.nanmax(mean['sp']))//5)*5
    levels_p = np.arange(min_p+10,max_p+5,5)
    
    CS = ax.contour(ds_sp.longitude,ds_sp.latitude,mean['sp'],colors='grey',zorder=0,levels = levels_p,linewidths=0.5,linestyles='--')   
    cs = ax.clabel(CS, fontsize=15, inline=1)
    
    ax.add_feature(cfeature.LAND, zorder=1, edgecolor='k',color='w')
    ax.add_feature(cfeature.COASTLINE,zorder=2,linewidth=0.8)
    ax.add_feature(cfeature.BORDERS,zorder=2,linewidth=0.2,alpha=0.5)
    ax.add_feature(cfeature.RIVERS,zorder=2,linewidth=3,alpha=0.5)
    ax.add_feature(cfeature.STATES,zorder=2,facecolor='w')
    ax.add_feature(cfeature.LAKES,zorder=2,linewidth=0.8,edgecolor='k',alpha=0.5)
         
    #Evaporation
    level_max = (np.round(np.nanmax(e_track))//step)*step
    level_min = 0
    levels_cb = np.arange(0,level_max+0.1,0.1)
    levels_ticks = np.arange(0,level_max+step,step)
        
    cb = ax.contourf(mean.coords['longitude'],mean.coords['latitude'],e_track,cmap=cmap_own, levels = levels_cb,zorder=3)
    cbar = plt.colorbar(cb,ticks=levels_ticks)
    cbar.ax.set_title('(mm d$^{-1}$)', fontsize=fontsizes,pad=pads)
    cbar.ax.tick_params(labelsize=fontsizes)
    
    #Add moisture flux   
    longitude_quiver = mean.coords['longitude'].values
    converted_longitudes = (longitude_quiver + 180) % 360 - 180
    quiver = ax.quiver(converted_longitudes,mean.coords['latitude'],mean_subset['mean_x'],mean_subset['mean_y'],zorder=4,transform=crs.PlateCarree(),regrid_shape=20,label='Moisture Flux')            
    qk = ax.quiverkey(quiver, 0.1, 1.05, 200, label='Moisture Flux', labelpos='N', labelcolor='k')
    qk.text.set_fontsize(fontsizes)
    
    #Add border of the configfile
    region.plot.contour(ax=ax, levels = [0.1,1], linewidths=3,colors='k',zorder=4)
    
    #Plot limit  
    index = np.where(region.values==1)
    latitude_sel = latitude[index[0]]
    longitude_sel = longitude[index[1]]
      
    ax.set_ylim(np.nanmin(latitude_sel)-50,np.nanmax(latitude_sel)+50)
    
    converted_longitudes = (longitude_sel + 180) % 360 - 180
    ax.set_xlim(np.nanmin(converted_longitudes)-50,np.nanmax(converted_longitudes)+50)
        
    #Grid
    gl = ax.gridlines(draw_labels=True,alpha=0.5, linestyle=':', color='k')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': fontsizes}
    gl.ylabel_style = {'size': fontsizes}
    
    # Save
    out_dir = Path(config.output_folder) / "figures"
    out_dir.mkdir(exist_ok=True, parents=True)
    fig.savefig(outdir+"/Evaporation_all_"+filename+".png", bbox_inches='tight', dpi=200)    
    

def moisture_fluxes(config_file,filename,data_path,outdir,step, fontsizes, pads):
    """
    
    A figure showing the moisture fluxes and their magnitude calculated by WAM-2Layers for the top and bottom layer
    
    Args:
      config_file (str): path to the .yaml configuration file to run the model. 
      filename (str): string to specify the output file.
      data_path (str): path to the location where the ERA5 (or other) data is stored. 
      Make sure to have data per year and per month. 
      
      outdir (str): path to an output directory to store the output config file. 
      step (int): stepsize of the ticks for the colourbar. 
      fontsizes (int): fontsizes. 
      pads (int): integer indicating the spacing between the labels and the axes. 
      
    Returns:
      two figures (.png) showing the mean moisture fluxes calculated by WAM-2Layers for the upper and lower layer of the model. 
        
    
    """

    #Load data
    config = Config.from_yaml(config_file)
    region = load_tagging_region(config)
    a_gridcell, lx, ly = get_grid_info(region)
    output_files = f"{config.output_folder}/*.nc"
    ds = xr.open_mfdataset(output_files, combine='nested', concat_dim='time')
    
    start = config.tracking_start_date
    end = config.tracking_end_date
            
    #Background data - isobars, wind field etc. 
    
    #Start date
    start_date = str(start)
    start_date = datetime.strptime(start_date, '%Y-%m-%d %H:%M:%S')
    
    #End date 
    end_date = str(end)
    end_date = datetime.strptime(end_date, '%Y-%m-%d %H:%M:%S')

            
    date_time = []
    one_month_delta = relativedelta(months=1)
    current_date = start_date
    while current_date <= end_date:
        string = current_date.strftime('%Y%m%d')
        date_time.append(datetime.strptime(string,'%Y%m%d'))
        current_date += one_month_delta #timedelta(month=1)
    
    months=[]
    years=[]
    for item in date_time:
        months.append(str(item)[5:7])
        years.append(str(item)[:4])
    
    
    #Surface pressure
    print('now loading p')
    
    all_sp = []
    for year in years:
        for month in months:
            ERA5_data_path = data_path+year+'/'+month+'/'
            subs = '*sp.nc'
            ds = xr.open_mfdataset(ERA5_data_path+subs)
            ds = ds / 100
            all_sp.append(ds)

    ds_sp = xr.concat(all_sp, dim="time")
    mean = ds_sp.mean(dim='time')
    
    print('p loaded')
    
    #Moisture fluxes
    path_input = f"{config.preprocessed_data_folder}/*.nc"
    ds = xr.open_mfdataset(path_input)
    
    start = config.tracking_start_date
    end = config.tracking_end_date
    subset = ds.sel(time=slice(start, end))
      
    mean_subset = subset.mean(dim='time')        
    mean_subset = mean_subset.assign(mean_lower=(np.sqrt((mean_subset['fy_lower'])**2 + (mean_subset['fx_lower'])**2)))
    mean_subset = mean_subset.assign(mean_upper=(np.sqrt((mean_subset['fy_upper'])**2 + (mean_subset['fx_upper'])**2)))
          
    #Get limits
    longitude = np.array(region.longitude.values)
    latitude = np.array(region.latitude.values)
    
    #own colourmap
    cmap_name = 'my_list'
    colors = ['w','#00008B','#800080','#E50000','#FFD700']  
    bin = 200
    cmap_own = LinearSegmentedColormap.from_list(cmap_name, colors, N=bin)

    # Make figure - lower
    fig = plt.figure(figsize=(16, 10))
    ax = fig.add_subplot(111, projection=crs.PlateCarree())
    
    #pressure
    min_p = (np.round(np.nanmean(mean['sp']))//5)*5
    max_p = (np.round(np.nanmax(mean['sp']))//5)*5
    levels_p = np.arange(min_p+10,max_p+5,5)
    
    CS = xr.plot.contour(mean['sp'],ax=ax,colors='grey',zorder=0,levels = levels_p,linewidths=0.5,linestyles='--')
    cs = ax.clabel(CS, fontsize=15, inline=1)
    
    #map features
    ax.add_feature(cfeature.LAND, zorder=1, edgecolor='k',color='w')
    ax.add_feature(cfeature.COASTLINE,zorder=2,linewidth=2)
    ax.add_feature(cfeature.BORDERS,zorder=2,linewidth=0.2,alpha=0.5)
    ax.add_feature(cfeature.RIVERS,zorder=2,linewidth=3,alpha=0.5)
    ax.add_feature(cfeature.STATES,zorder=2,facecolor='w')
    ax.add_feature(cfeature.LAKES,zorder=2,linewidth=0.8,edgecolor='k',alpha=0.5)
      
    #Moisture flux 
    level_max = (np.round(np.nanmax(mean_subset.mean_lower))//step)*step
    level_min = 0
    levels_cb = np.arange(0,level_max+1,1)
    levels_ticks = np.arange(0,level_max+step,step)
    
    cb = ax.contourf(mean.coords['longitude'],mean.coords['latitude'],mean_subset.mean_lower,cmap=cmap_own, levels = levels_cb,zorder=3,alpha=0.4)
    cbar = plt.colorbar(cb,ticks=levels_ticks)
    cbar.ax.set_title('(kg m$^{-1}$ s$^{-1}$)', fontsize=fontsizes,pad=pads)
    cbar.ax.tick_params(labelsize=fontsizes)
    
    #Add arrows
    longitude_quiver = mean.coords['longitude'].values
    converted_longitudes = (longitude_quiver + 180) % 360 - 180
    quiver = ax.quiver(converted_longitudes,mean_subset.coords['latitude'],mean_subset['fx_lower'],mean_subset['fy_lower'],transform=crs.PlateCarree(),regrid_shape=20,label='Moisture Flux',zorder=3)         
    
    qk = ax.quiverkey(quiver, 0.1, 1.05, 200, label='Moisture Flux', labelpos='N', labelcolor='k')
    qk.text.set_fontsize(fontsizes)
    
    #Add border of the configfile
    region.plot.contour(ax=ax, levels = [0.1,1], linewidths=3,colors='k',zorder=4)
    
    #Plot limit
    index = np.where(region.values==1)
    latitude_sel = latitude[index[0]]
    longitude_sel = longitude[index[1]]
      
    ax.set_ylim(np.nanmin(latitude_sel)-50,np.nanmax(latitude_sel)+50)
    
    converted_longitudes = (longitude_sel + 180) % 360 - 180
    ax.set_xlim(np.nanmin(converted_longitudes)-50,np.nanmax(converted_longitudes)+50)
    
    #Grid
    gl = ax.gridlines(draw_labels=True,alpha=0.5, linestyle=':', color='k')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': fontsizes}
    gl.ylabel_style = {'size': fontsizes}
    
    # Save
    out_dir = Path(config.output_folder) / "figures"
    out_dir.mkdir(exist_ok=True, parents=True)
    fig.savefig(outdir+"/Moisture_flux_lower_"+filename+".png", bbox_inches='tight', dpi=200)
    
    # Make figure - upper
    fig = plt.figure(figsize=(16, 10))
    ax = fig.add_subplot(111, projection=crs.PlateCarree())
    
    #pressure
    min_p = (np.round(np.nanmean(mean['sp']))//5)*5
    max_p = (np.round(np.nanmax(mean['sp']))//5)*5
    levels_p = np.arange(min_p+10,max_p+5,5)
    
    CS = xr.plot.contour(mean['sp'],ax=ax,colors='grey',zorder=0,levels = levels_p,linewidths=0.5,linestyles='--')
    cs = ax.clabel(CS, fontsize=15, inline=1)
    
    #map features
    ax.add_feature(cfeature.LAND, zorder=1, edgecolor='k',color='w')
    ax.add_feature(cfeature.COASTLINE,zorder=2,linewidth=2)
    ax.add_feature(cfeature.BORDERS,zorder=2,linewidth=0.2,alpha=0.5)
    ax.add_feature(cfeature.RIVERS,zorder=2,linewidth=3,alpha=0.5)
    ax.add_feature(cfeature.STATES,zorder=2,facecolor='w')
    ax.add_feature(cfeature.LAKES,zorder=2,linewidth=0.8,edgecolor='k',alpha=0.5)
    
        
    #Moisture flux 
    level_max = (np.round(np.nanmax(mean_subset.mean_lower))//step)*step
    level_min = 0
    levels_cb = np.arange(0,level_max+1,1)
    levels_ticks = np.arange(0,level_max+step,step)
    
    cb = ax.contourf(mean.coords['longitude'],mean.coords['latitude'],mean_subset.mean_upper,cmap=cmap_own, levels = levels_cb,zorder=3,alpha=0.4)
    cbar = plt.colorbar(cb,ticks=levels_ticks)
    cbar.ax.set_title('(kg m$^{-1}$ s$^{-1}$)', fontsize=fontsizes,pad=pads)
    cbar.ax.tick_params(labelsize=fontsizes)
    
    #Add arrows
    longitude_quiver = mean.coords['longitude'].values
    converted_longitudes = (longitude_quiver + 180) % 360 - 180
    quiver = ax.quiver(converted_longitudes,mean_subset.coords['latitude'],mean_subset['fx_upper'],mean_subset['fy_upper'],zorder=4,transform=crs.PlateCarree(),regrid_shape=20,label='Moisture Flux')         
    
    qk = ax.quiverkey(quiver, 0.1, 1.05, 200, label='Moisture Flux', labelpos='N', labelcolor='k')
    qk.text.set_fontsize(fontsizes)
    
    #Add border of the configfile
    region.plot.contour(ax=ax, levels = [1], linewidths=3,colors='k',zorder=4)
    
    #Plot limit
    index = np.where(region.values==1)
    latitude_sel = latitude[index[0]]
    longitude_sel = longitude[index[1]]
      
    ax.set_ylim(np.nanmin(latitude_sel)-50,np.nanmax(latitude_sel)+50)
    
    converted_longitudes = (longitude_sel + 180) % 360 - 180
    ax.set_xlim(np.nanmin(converted_longitudes)-50,np.nanmax(converted_longitudes)+50)

        
    #Grid
    gl = ax.gridlines(draw_labels=True,alpha=0.5, linestyle=':', color='k')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': fontsizes}
    gl.ylabel_style = {'size': fontsizes}
    
    # Save
    out_dir = Path(config.output_folder) / "figures"
    out_dir.mkdir(exist_ok=True, parents=True)
    fig.savefig(outdir+"/Moisture_flux_upper_"+filename+".png", bbox_inches='tight', dpi=200)