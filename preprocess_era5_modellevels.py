import os

import pandas as pd
import xarray as xr
import yaml
import numpy as np

from preprocessing import (get_grid_info, get_stable_fluxes,
                           get_vertical_transport)


# Set constants
g = 9.80665  # [m/s2]
density_water = 1000  # [kg/m3]

# Select model levels
modellevels = np.append( np.arange(1,100,10), np.arange(101,138,1)) #[1,20,40,60,80,100,110,120,125,130,131,132,133,134,135,136,137]
print('Number of model levels:',len(modellevels))

# Read case configuration
with open("cases/era5_2013.yaml") as f:
    config = yaml.safe_load(f)

def load_surface_data(variable, date):
    """Load data for given variable and date."""
    filename = f"FloodCase_201305_{variable}.nc"
    filepath = os.path.join(config["input_folder"], filename)
    da = xr.open_dataset(filepath)[variable]

    # Include midnight of the next day (if available)
    extra = date + pd.Timedelta(days=1)
    return da.sel(time=slice(date, extra))

def load_modellevel_data(variable, date):
    """Load model level data for given variable and date."""
    filename = f"FloodCase_201305_ml_{variable}.nc"
    filepath = os.path.join(config["input_folder"], filename)
    da = xr.open_dataset(filepath)[variable]

    # Include midnight of the next day (if available)
    extra = date + pd.Timedelta(days=1)
    return da.sel(time=slice(date, extra)).sel(lev=modellevels)

datelist = pd.date_range(
start=config["preprocess_start_date"], end=config["preprocess_end_date"], freq="d", inclusive="left"
)

for date in datelist[:1]:
    print(date)

    # Load data
    u = load_modellevel_data("u", date)
    v = load_modellevel_data("v", date)
    q = load_modellevel_data("q", date)
    sp = load_surface_data("sp", date) # in Pa
    evap = load_surface_data("e", date)
    cp = load_surface_data("cp", date)
    lsp = load_surface_data("lsp", date)
    tcw = load_surface_data("tcw", date) # kg/m2
    tcrw = load_surface_data("tcrw", date) # kg/m2
    precip = cp + lsp

    # Get grid info
    lat = u.latitude.values
    lon = u.longitude.values
    a_gridcell, l_ew_gridcell, l_mid_gridcell = get_grid_info(lat, lon)

    # Calculate volumes
    evap *= a_gridcell  # m3
    precip *= a_gridcell  # m3

    # Convert model levels to pressure values
    filenamecsv = 'tableERA5model_to_pressure.csv'
    df = pd.read_csv(os.path.join(config["input_folder"],filenamecsv))
    a = df['a [Pa]'].to_xarray().rename(index='lev').sel(lev=modellevels)
    b = df['b'].to_xarray().rename(index='lev').sel(lev=modellevels)
    
    sp_modellevels2 = sp.expand_dims({'lev':modellevels},axis=1)
    p_modellevels = a + b*sp_modellevels2 # in Pa

    #calculate the difference between the pressure levels
    dp_modellevels = p_modellevels.diff(dim="lev") # in Pa

    # Determine the fluxes and states
    fa_e = u * q * dp_modellevels / g  # eastward atmospheric moisture flux
    fa_n = v * q * dp_modellevels / g  # northward atmospheric moisture flux
    cwv = q * dp_modellevels / g * a_gridcell[np.newaxis,np.newaxis,:] / density_water  # column water vapor (m3)

    # Split in 2 layers 
    boundary = 111 # set boundary model level 111
    
    # Integrate fluxes and states to upper and lower layer
    lower_layer = dp_modellevels.lev > boundary 
    upper_layer = ~lower_layer

    # Integrate fluxes and state
    fa_e_lower = fa_e.where(lower_layer).sum(dim="lev")
    fa_n_lower = fa_n.where(lower_layer).sum(dim="lev")
    w_lower = cwv.where(lower_layer).sum(dim="lev")

    fa_e_upper = fa_e.where(upper_layer).sum(dim="lev")
    fa_n_upper = fa_n.where(upper_layer).sum(dim="lev")
    w_upper = cwv.where(upper_layer).sum(dim="lev")
   
    print(
        "Check calculation water vapor over two layers, this value should be zero:",
        (cwv.sum(dim="lev") - (w_upper + w_lower)).sum().values,
    )

    tcwm3 = tcw* a_gridcell[np.newaxis,:] / density_water # m3
    
    print(
        "Check calculation water vapor with column water vapour, this value should be close to zero:",
        (cwv.sum(dim="lev") - (tcwm3)).sum().values,
    )
    # is niet hetzelfde.. hoe erg is dat?
    
    # Change units to m3, based on target frequency (not incoming frequency!)
    target_freq = config['target_frequency']
    total_seconds = pd.Timedelta(target_freq).total_seconds()
    fa_e_upper *= total_seconds * (l_ew_gridcell / density_water)
    fa_e_lower *= total_seconds * (l_ew_gridcell / density_water)
    fa_n_upper *= total_seconds * (l_mid_gridcell[None, :, None] / density_water)
    fa_n_lower *= total_seconds * (l_mid_gridcell[None, :, None] / density_water)

    # Put data on a smaller time step...
    time = w_upper.time.values
    newtime = pd.date_range(time[0], time[-1], freq="15Min")[:-1]
    w_upper = w_upper.interp(time=newtime).values
    w_lower = w_lower.interp(time=newtime).values

    # ... fluxes on the edges instead of midpoints
    newtime = newtime[:-1] + pd.Timedelta("6Min") / 2
    fa_e_upper = fa_e_upper.interp(time=newtime).values
    fa_n_upper = fa_n_upper.interp(time=newtime).values

    fa_e_lower = fa_e_lower.interp(time=newtime).values
    fa_n_lower = fa_n_lower.interp(time=newtime).values

    precip = (precip.reindex(time=newtime, method="bfill") / 4).values
    evap = (evap.reindex(time=newtime, method="bfill") / 4).values

    # Stabilize horizontal fluxes
    fa_e_upper, fa_n_upper = get_stable_fluxes(fa_e_upper, fa_n_upper, w_upper)
    fa_e_lower, fa_n_lower = get_stable_fluxes(fa_e_lower, fa_n_lower, w_lower)

    # Determine the vertical moisture flux
    fa_vert = get_vertical_transport(
        fa_e_upper,
        fa_e_lower,
        fa_n_upper,
        fa_n_lower,
        evap,
        precip,
        w_upper,
        w_lower,
    )

    # Save preprocessed data
    # Note: fluxes (dim: time) are at the edges of the timesteps,
    # while states (dim: time2) are at the midpoints and include next midnight
    # so the first state from day 2 will overlap with the last flux from day 1
    filename = f"{date.strftime('%Y-%m-%d')}_fluxes_storages.nc"
    output_path = os.path.join(config["preprocessed_data_folder"], filename)
    xr.Dataset(
        {  # TODO: would be nice to add coordinates and units as well
            "fa_e_upper": (["time", "lat", "lon"], fa_e_upper),
            "fa_n_upper": (["time", "lat", "lon"], fa_n_upper),
            "fa_e_lower": (["time", "lat", "lon"], fa_e_lower),
            "fa_n_lower": (["time", "lat", "lon"], fa_n_lower),
            "w_upper": (["time2", "lat", "lon"], w_upper),
            "w_lower": (["time2", "lat", "lon"], w_lower),
            "fa_vert": (["time", "lat", "lon"], fa_vert),
            "evap": (["time", "lat", "lon"], evap),
            "precip": (["time", "lat", "lon"], precip),
        }
    ).to_netcdf(output_path)
