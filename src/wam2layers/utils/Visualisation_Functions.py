# coding: cp1252
''' 

Visualisation Functions WAM-2Layers

Version 1.0 | 25-06-2023

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
from wam2layers.preprocessing.shared import get_grid_info
from wam2layers.tracking.backtrack import (input_path, load_region,
                                           output_path, parse_config)
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
    config = parse_config(config_file)
    region = load_region(config)
    
    #Load precipitation
    input_files = f"{config['preprocessed_data_folder']}/*.nc"
    ds_p = xr.open_mfdataset(input_files, combine='nested', concat_dim='time')
    start = config["event_start_date"]
    end = config["event_end_date"]
    subset = ds_p.precip.sel(time=slice(start, end))
    precip = (subset * region * 3600).sum('time').compute()
    
    
    #Load tracked moisture
    a_gridcell, lx, ly = get_grid_info(region)
    output_files = f"{config['output_folder']}/*.nc"
    ds_e = xr.open_mfdataset(output_files, combine='nested', concat_dim='time')
    e_track = ds_e.e_track.sum('time').compute() * 1000 / a_gridcell [:, None]
        
    #Get moisture flux data
    path_input = f"{config['preprocessed_data_folder']}/*.nc"
    ds = xr.open_mfdataset(path_input)
    
    start = config["track_start_date"]
    end = config["track_end_date"]
    subset = ds.sel(time=slice(start, end))
    
    mean_subset = ds.sum(dim='time')
       
    
    #calculate evaporation and precipitation ratio
       
    mean_subset['evap_masked'] = mean_subset['evap'].where(region, drop=True) * 1000
    
    mean_subset['precip_masked'] = mean_subset['precip'].where(region, drop=True) * 1000
    
    mean_subset['tracked_masked'] = e_track.where(region, drop=True)
    
    mean_subset = mean_subset.assign(e_r=(mean_subset['tracked_masked'].sum() / mean_subset['evap_masked'].sum())*100)
    mean_subset = mean_subset.assign(p_r=(mean_subset['tracked_masked'].sum() / mean_subset['precip_masked'].sum())*100)  
  
    precip_sum = precip.sum()
    e_track_sum = e_track.sum()
    
    diff = abs(precip_sum - e_track_sum)
    
    
    #Write to file
    logfile = open(outdir+"/Logfile_"+filename+".txt",'w')
    logfile.write('precip_sum: '+str(float(precip_sum))+'\n')
    logfile.write('e_track_sum: '+str(float(e_track_sum))+'\n')
    logfile.write('abs_diff: '+str(float(diff))+'\n')
    logfile.write('e_r: '+str(float(mean_subset['e_r']))+'\n')
    logfile.write('p_r: '+str(float(mean_subset['p_r']))+'\n')
    
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
    config = parse_config(config_file)
    region = load_region(config)
    
    input_files = f"{config['preprocessed_data_folder']}/*.nc"
    ds = xr.open_mfdataset(input_files, combine='nested', concat_dim='time')
    start = config["event_start_date"]
    end = config["event_end_date"]
        
    #Load precipitaiton
    subset = ds.precip.sel(time=slice(start, end))
    precip = (subset * region * 3600).sum('time').compute()
    precip = np.where(precip>1,precip,np.nan)
    
    #Get limits
    index = np.where(region==1)
    longitude = np.array(region.longitude)
    latitude = np.array(region.latitude)
    longitude = longitude[index[1]]
    latitude = latitude[index[0]]

           
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
    cbar.ax.set_title('[mm]', fontsize=fontsizes,pad=pads)
    cbar.ax.tick_params(labelsize=fontsizes)    
    
    #Add border of the configfile
    region.plot.contour(ax=ax, levels = [1], linewidths=3,colors='k',zorder=4)
    
    #Plot limit
    x_min = round(np.nanmin(longitude)/10)*10
    x_max = round(np.nanmax(longitude)/10)*10
    
    y_min = round(np.nanmin(latitude)/10)*10
    y_max = round(np.nanmax(latitude)/10)*10
    
    
    ax.set_xlim(x_min-10,x_max+10)
    ax.set_ylim(y_min-10,y_max+10)
        
    #Add minimap
    cmap_name = 'my_list'
    colors = ['tab:blue','tab:blue']  
    bin = 50
    cmap_own = LinearSegmentedColormap.from_list(cmap_name, colors, N=bin)
    
    
    ax_mini = fig.add_axes([0.15, 0.15, 0.2, 0.2], projection=crs.PlateCarree(),facecolor='None') #x,y,w,h
    ax_mini.set_extent([x_min-50,x_max+50, y_min-50,y_max+50], crs=crs.PlateCarree())
    ax_mini.coastlines(alpha=0.5)
    ax_mini.add_feature(cfeature.COASTLINE,alpha=0.5)
    ax_mini.add_feature(cfeature.BORDERS,alpha=0.5)
    region.plot.contour(ax=ax_mini, levels = [1], linewidths=3,colors='k')
    
    #Grid
    gl = ax.gridlines(draw_labels=True,alpha=0.5, linestyle=':', color='k')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlocator = plt.FixedLocator(np.arange(x_min-10,x_max+20,10))
    gl.ylocator = plt.FixedLocator(np.arange(y_min-10,y_max+20,10))
    gl.xlabel_style = {'size': fontsizes}
    gl.ylabel_style = {'size': fontsizes}
    
    # Export figure
    out_dir = Path(config["output_folder"]) / "figures"
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
    config = parse_config(config_file)
    region = load_region(config)
    a_gridcell, lx, ly = get_grid_info(region)
    output_files = f"{config['output_folder']}/*.nc"
    ds = xr.open_mfdataset(output_files, combine='nested', concat_dim='time')
    
    start = config["track_start_date"]
    end = config["track_end_date"]
    
    #Load evaporation
    e_track = ds.e_track.sum('time').compute() * 1000 / a_gridcell [:, None]
    e_track = np.where(e_track>1,e_track,np.nan)
        
    #Background data - isobars, wind field etc. 

    
    #Start date
    start_date = start
    start_date = datetime.strptime(start_date, '%Y%m%d')
    
    #End date 
    end_date = end
    end_date = datetime.strptime(end_date, '%Y%m%d')

            
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
    
    #v
    print('now loading v')
    all_v = []
    for year in years:
      for month in months:
        ERA5_data_path = data_path+year+'/'+month+'/'
        subs = '*ml_v.nc'
        ds = xr.open_mfdataset(ERA5_data_path+subs)
        level_v = ds.sel(level=115)
        all_v.append(level_v)
      
    ds_v = xr.concat(all_v, dim="time")
    mean_v = ds_v.mean(dim='time')
    
    print('v loaded')
    
    #u
    print('now loading u')
    all_u = []
    for year in years:
      for month in months:
        ERA5_data_path = data_path+year+'/'+month+'/'
        subs = '*ml_u.nc'
        ds = xr.open_mfdataset(ERA5_data_path+subs)
        level_u = ds.sel(level=115)
        all_u.append(level_u)
      
    ds_u = xr.concat(all_u, dim="time")
    mean_u = ds_u.mean(dim='time')
    
    print('u loaded')
    
    #Moisture fluxes
    path_input = f"{config['preprocessed_data_folder']}/*.nc"
    ds = xr.open_mfdataset(path_input)
    
    start = config["track_start_date"]
    end = config["track_end_date"]
    subset = ds.sel(time=slice(start, end))
      
    mean_subset = subset.mean(dim='time')
        
    mean_subset = mean_subset.assign(mean_x=(mean_subset['fx_upper'] + mean_subset['fx_lower'])/2)
    
    mean_subset = mean_subset.assign(mean_y=(mean_subset['fy_upper'] + mean_subset['fy_lower'])/2)
    
    #Get limits
    index = np.where(region==1)
    longitude = np.array(region.longitude)
    latitude = np.array(region.latitude)
    longitude = longitude[index[1]]
    latitude = latitude[index[0]]
    
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
    
    CS = xr.plot.contour(mean['sp'],ax=ax,colors='grey',zorder=0,levels = levels_p,linewidths=0.5,linestyles='--')   
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
    levels_cb = np.arange(0,level_max+1,1)
    levels_ticks = np.arange(0,level_max+step,step)
    
    cb = ax.contourf(mean.coords['longitude'],mean.coords['latitude'],e_track,cmap=cmap_own, levels = levels_cb,zorder=3)
    cbar = plt.colorbar(cb,ticks=levels_ticks)
    cbar.ax.set_title('[mm d$^{-1}$]', fontsize=fontsizes,pad=pads)
    cbar.ax.tick_params(labelsize=fontsizes)
    
    #Add moisture flux
    quiver = ax.quiver(mean_subset.coords['longitude'],mean_subset.coords['latitude'],mean_subset['mean_x'],mean_subset['mean_y'],zorder=4,transform=crs.PlateCarree(),regrid_shape=20,label='Moisture Flux')         
    
    qk = ax.quiverkey(quiver, 0.1, 1.05, 200, label='Moisture Flux', labelpos='N', labelcolor='k')
    qk.text.set_fontsize(fontsizes)
    
    
    #Add border of the configfile
    region.plot.contour(ax=ax, levels = [1], linewidths=3,colors='k',zorder=4)
    
    #Plot limit
    x_min = round(np.nanmin(longitude)/10)*10
    x_max = round(np.nanmax(longitude)/10)*10
    
    y_min = round(np.nanmin(latitude)/10)*10
    y_max = round(np.nanmax(latitude)/10)*10
    
    ax.set_xlim(x_min-50,x_max+50)
    ax.set_ylim(y_min-50,y_max+50)
        
    #Grid
    gl = ax.gridlines(draw_labels=True,alpha=0.5, linestyle=':', color='k')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlocator = plt.FixedLocator(np.arange(x_min-50,x_max+70,20))
    gl.ylocator = plt.FixedLocator(np.arange(y_min-50,y_max+70,20))
    gl.xlabel_style = {'size': fontsizes}
    gl.ylabel_style = {'size': fontsizes}
    
    # Save
    out_dir = Path(config["output_folder"]) / "figures"
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
    config = parse_config(config_file)
    region = load_region(config)
    a_gridcell, lx, ly = get_grid_info(region)
    output_files = f"{config['output_folder']}/*.nc"
    ds = xr.open_mfdataset(output_files, combine='nested', concat_dim='time')
    
    start = config["track_start_date"]
    end = config["track_end_date"]
            
    #Background data - isobars, wind field etc. 
    
    #Start date
    start_date = start
    start_date = datetime.strptime(start_date, '%Y%m%d')
    
    #End date 
    end_date = end
    end_date = datetime.strptime(end_date, '%Y%m%d')

            
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
    
    #v
    print('now loading v')
    all_v = []
    for year in years:
      for month in months:
        ERA5_data_path = data_path+year+'/'+month+'/'
        subs = '*ml_v.nc'
        ds = xr.open_mfdataset(ERA5_data_path+subs)
        level_v = ds.sel(level=115)
        all_v.append(level_v)
      
    ds_v = xr.concat(all_v, dim="time")
    mean_v = ds_v.mean(dim='time')
    
    print('v loaded')
    
    #u
    print('now loading u')
    all_u = []
    for year in years:
      for month in months:
        ERA5_data_path = data_path+year+'/'+month+'/'
        subs = '*ml_u.nc'
        ds = xr.open_mfdataset(ERA5_data_path+subs)
        level_u = ds.sel(level=115)
        all_u.append(level_u)
      
    ds_u = xr.concat(all_u, dim="time")
    mean_u = ds_u.mean(dim='time')
    
    print('u loaded')
    
    #Moisture fluxes
    path_input = f"{config['preprocessed_data_folder']}/*.nc"
    ds = xr.open_mfdataset(path_input)
    
    start = config["track_start_date"]
    end = config["track_end_date"]
    subset = ds.sel(time=slice(start, end))
      
    mean_subset = subset.mean(dim='time')
        
    mean_subset = mean_subset.assign(mean_lower=(np.sqrt((mean_subset['fy_lower'])**2 + (mean_subset['fx_lower'])**2)))
    
    mean_subset = mean_subset.assign(mean_upper=(np.sqrt((mean_subset['fy_upper'])**2 + (mean_subset['fx_upper'])**2)))
       
    
    #Get limits
    index = np.where(region==1)
    longitude = np.array(region.longitude)
    latitude = np.array(region.latitude)
    longitude = longitude[index[1]]
    latitude = latitude[index[0]]
    
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
    cbar.ax.set_title('[kg m$^{-1}$ s$^{-1}$]', fontsize=fontsizes,pad=pads)
    cbar.ax.tick_params(labelsize=fontsizes)
    
    #Add arrows
    quiver = ax.quiver(mean_subset.coords['longitude'],mean_subset.coords['latitude'],mean_subset['fx_lower'],mean_subset['fy_lower'],zorder=4,transform=crs.PlateCarree(),regrid_shape=20,label='Moisture Flux')         
    
    qk = ax.quiverkey(quiver, 0.1, 1.05, 200, label='Moisture Flux', labelpos='N', labelcolor='k')
    qk.text.set_fontsize(fontsizes)
    
    #Add border of the configfile
    region.plot.contour(ax=ax, levels = [1], linewidths=3,colors='k',zorder=4)
    
    #Plot limit
    x_min = round(np.nanmin(longitude)/10)*10
    x_max = round(np.nanmax(longitude)/10)*10
    
    y_min = round(np.nanmin(latitude)/10)*10
    y_max = round(np.nanmax(latitude)/10)*10
    
    ax.set_xlim(x_min-50,x_max+50)
    ax.set_ylim(y_min-50,y_max+50)
        
    #Grid
    gl = ax.gridlines(draw_labels=True,alpha=0.5, linestyle=':', color='k')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlocator = plt.FixedLocator(np.arange(x_min-50,x_max+70,20))
    gl.ylocator = plt.FixedLocator(np.arange(y_min-50,y_max+70,20))
    gl.xlabel_style = {'size': fontsizes}
    gl.ylabel_style = {'size': fontsizes}
    
    # Save
    out_dir = Path(config["output_folder"]) / "figures"
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
    cbar.ax.set_title('[kg m$^{-1}$ s$^{-1}$]', fontsize=fontsizes,pad=pads)
    cbar.ax.tick_params(labelsize=fontsizes)
    
    #Add arrows
    quiver = ax.quiver(mean_subset.coords['longitude'],mean_subset.coords['latitude'],mean_subset['fx_upper'],mean_subset['fy_upper'],zorder=4,transform=crs.PlateCarree(),regrid_shape=20,label='Moisture Flux')         
    
    qk = ax.quiverkey(quiver, 0.1, 1.05, 200, label='Moisture Flux', labelpos='N', labelcolor='k')
    qk.text.set_fontsize(fontsizes)
    
    #Add border of the configfile
    region.plot.contour(ax=ax, levels = [1], linewidths=3,colors='k',zorder=4)
    
    #Plot limit
    x_min = round(np.nanmin(longitude)/10)*10
    x_max = round(np.nanmax(longitude)/10)*10
    
    y_min = round(np.nanmin(latitude)/10)*10
    y_max = round(np.nanmax(latitude)/10)*10
    
    ax.set_xlim(x_min-50,x_max+50)
    ax.set_ylim(y_min-50,y_max+50)
        
    #Grid
    gl = ax.gridlines(draw_labels=True,alpha=0.5, linestyle=':', color='k')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlocator = plt.FixedLocator(np.arange(x_min-50,x_max+70,20))
    gl.ylocator = plt.FixedLocator(np.arange(y_min-50,y_max+70,20))
    gl.xlabel_style = {'size': fontsizes}
    gl.ylabel_style = {'size': fontsizes}
    
    # Save
    out_dir = Path(config["output_folder"]) / "figures"
    out_dir.mkdir(exist_ok=True, parents=True)
    fig.savefig(outdir+"/Moisture_flux_upper_"+filename+".png", bbox_inches='tight', dpi=200)