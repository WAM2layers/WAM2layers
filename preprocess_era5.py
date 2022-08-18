import numpy as np
import pandas as pd
import xarray as xr
import yaml
from pathlib import Path

from preprocessing import get_grid_info, get_stable_fluxes, get_vertical_transport


# Set constants
g = 9.80665  # [m/s2]
density_water = 1000  # [kg/m3]

# Read case configuration
with open("cases/era5_2013.yaml") as f:
    config = yaml.safe_load(f)

# Create the preprocessed data folder if it does not exist yet
output_dir = Path(config["preprocessed_data_folder"]).expanduser()
output_dir.mkdir(exist_ok=True, parents=True)


def load_data(variable, date):
    """Load data for given variable and date."""
    filepath = Path(config["input_folder"]) / f"FloodCase_201305_{variable}.nc"
    da = xr.open_dataset(filepath)[variable]

    # Include midnight of the next day (if available)
    extra = date + pd.Timedelta(days=1)
    return da.sel(time=slice(date, extra))


datelist = pd.date_range(
    start=config["preprocess_start_date"],
    end=config["preprocess_end_date"],
    freq="d",
    inclusive="left",
)

def calc_q2m(d2m,ps):
    '''
    Calculate the specific humidity from (surface) pressure and
    dew point temperature

    See further details at eq. 7.4 and 7.5 (Page 102) of:
    https://www.ecmwf.int/en/elibrary/20198-ifs-documentation-cy47r3-part-iv-physical-processes
    '''
    R_d =287.0597
    R_v=461.5250
    a1=611.21
    a3=17.502
    a4=32.19
    T0=273.15

    #Calculation of E saturation water vapour pressure from Teten's formula
    E=a1*np.exp(a3*(d2m-T0)/(d2m-a4))

    #Specific humidity
    q2m=(R_d/R_v)*E/(ps-((1-R_d/R_v)*E))

    return q2m

for date in datelist[:]:
    print(date)

    # Load data
    q = load_data("q", date)  # in kg kg-1
    u = load_data("u", date)  # in m/s
    v = load_data("v", date)  # in m/s
    sp = load_data("sp", date)  # in Pa
    evap = load_data("e", date)  # in m (accumulated hourly)
    cp = load_data("cp", date)  # convective precipitation in m (accumulated hourly)
    lsp = load_data("lsp", date)  # large scale precipitation in m (accumulated hourly)
    precip = cp + lsp
    tcw = load_data("tcw", date) # kg/m2
    d2m = load_data("d2m", date) # Dew point in K
    q_surf = calc_q2m(d2m,sp) # kg kg-1
    u_surf = load_data("u10",date) # in m/s
    v_surf = load_data("v10",date) # in m/s

    # Get grid info
    time = u.time.values
    lat = u.latitude.values
    lon = u.longitude.values
    a_gridcell, l_ew_gridcell, l_mid_gridcell = get_grid_info(lat, lon)

    # Calculate volumes
    evap *= a_gridcell  # m3
    precip *= a_gridcell  # m3

    # Transfer negative (originally positive) values of evap to precip
    precip = np.maximum(precip, 0) + np.maximum(evap, 0)
    # Change sign convention to all positive,
    evap = np.abs(np.minimum(evap, 0))

    # Create full pressure array (including top of atmosphere and surface pressure, 1100 is a dummy value)
    levels = np.append(0, np.append(np.array(q.level), [1100, 1100])) * 100 # Pa
    levels_reshaped = np.reshape(levels, [1, len(levels), 1, 1])
    p_temp = np.tile(levels_reshaped, (len(time), 1, len(lat), len(lon)))
    sp_reshaped = np.reshape(np.array(sp), [len(time), 1, len(lat), len(lon)])
    sp_temp = np.tile(sp_reshaped, (1, len(levels), 1, 1))

    # find the location of the surface
    above_surface = p_temp < sp_temp
    surf_loc = np.sum(above_surface, axis=1)

    # find the location of the boundary
    p_boundary = 0.72878581 * np.array(sp) + 7438.803223
    p_boundary_reshaped = np.reshape(np.array(p_boundary), [len(time), 1, len(lat), len(lon)])
    p_boundary_temp = np.tile(p_boundary_reshaped, (1, len(levels), 1, 1))
    above_boundary = (p_temp < p_boundary_temp)
    boundary_loc = np.sum(above_boundary, axis=1)

    # Create pressure fields: [time, top-upper_levels-boundary-lower_levels-surface:, latitude, longitude]
    # TODO: the loop below works, but is super slow, check if this loop is necessary (if removed possibly keep stored for debugging)
    p_full = np.zeros(p_temp.shape)
    q_full = np.zeros(p_temp.shape) # assume for p=0, q=0
    u_full = np.zeros(p_temp.shape)
    u_full[:,0,:,:] = u[:,0,:,:] # assume for p=0, u is equal to u at highest level in the atmosphere
    v_full = np.zeros(p_temp.shape)
    v_full[:,0,:,:] = v[:,0,:,:] # assume for p=0, u is equal to u at highest level in the atmosphere
    for t in range(len(time)):
        for i in range(len(lat)):
            for j in range(len(lon)):
                p_full[t,1:boundary_loc[t,i,j],i,j] = p_temp[t,1:boundary_loc[t,i,j],i,j]
                p_full[t,boundary_loc[t,i,j],i,j] = p_boundary[t,i,j]
                p_full[t,boundary_loc[t,i,j]+1:surf_loc[t,i,j]+1,i,j] = p_temp[t,boundary_loc[t,i,j]:surf_loc[t,i,j],i,j]
                p_full[t,surf_loc[t,i,j]+1:,i,j] = sp[t,i,j]

                q_full[t,1:boundary_loc[t,i,j],i,j] = q[t,0:boundary_loc[t,i,j]-1,i,j]
                q_full[t,boundary_loc[t,i,j],i,j] = q[t,boundary_loc[t,i,j]-2,i,j] + \
                                                    ( (q[t,boundary_loc[t,i,j]-1,i,j] - q[t,boundary_loc[t,i,j]-2,i,j]) \
                                                     * (p_boundary[t,i,j] - p_full[t,boundary_loc[t,i,j]-1,i,j]) ) \
                                                        / (p_full[t,boundary_loc[t,i,j]+1,i,j] - p_full[t,boundary_loc[t,i,j]-1,i,j])
                q_full[t,boundary_loc[t,i,j]+1:surf_loc[t,i,j]+1,i,j] = q[t,boundary_loc[t,i,j]-1:surf_loc[t,i,j]-1,i,j]
                q_full[t,surf_loc[t,i,j]:,i,j] = q_surf[t,i,j]

                u_full[t,1:boundary_loc[t,i,j],i,j] = u[t,0:boundary_loc[t,i,j]-1,i,j]
                u_full[t,boundary_loc[t,i,j],i,j] = u[t,boundary_loc[t,i,j]-2,i,j] + \
                                                    ( (u[t,boundary_loc[t,i,j]-1,i,j] - u[t,boundary_loc[t,i,j]-2,i,j]) \
                                                     * (p_boundary[t,i,j] - p_full[t,boundary_loc[t,i,j]-1,i,j]) ) \
                                                        / (p_full[t,boundary_loc[t,i,j]+1,i,j] - p_full[t,boundary_loc[t,i,j]-1,i,j])
                u_full[t,boundary_loc[t,i,j]+1:surf_loc[t,i,j]+1,i,j] = u[t,boundary_loc[t,i,j]-1:surf_loc[t,i,j]-1,i,j]
                u_full[t,surf_loc[t,i,j]:,i,j] = u_surf[t,i,j]

                v_full[t,1:boundary_loc[t,i,j],i,j] = v[t,0:boundary_loc[t,i,j]-1,i,j]
                v_full[t,boundary_loc[t,i,j],i,j] = v[t,boundary_loc[t,i,j]-2,i,j] + \
                                                    ( (v[t,boundary_loc[t,i,j]-1,i,j] - v[t,boundary_loc[t,i,j]-2,i,j]) \
                                                     * (p_boundary[t,i,j] - p_full[t,boundary_loc[t,i,j]-1,i,j]) ) \
                                                        / (p_full[t,boundary_loc[t,i,j]+1,i,j] - p_full[t,boundary_loc[t,i,j]-1,i,j])
                v_full[t,boundary_loc[t,i,j]+1:surf_loc[t,i,j]+1,i,j] = v[t,boundary_loc[t,i,j]-1:surf_loc[t,i,j]-1,i,j]
                v_full[t,surf_loc[t,i,j]:,i,j] = v_surf[t,i,j]

                # TODO, if we keep this code like this, create a function for the manipulation of q, u, and v (uses the same code)

    # TODO: add a quick check whether the new pressure level list is strictly increasing, otherwise give an error!

    # Interpolate to new levels
    midpoints = 0.5 * (p_full[:,:-1,:,:] + p_full[:,1:,:,:]) # TODO: not used, but could be to used when converted to x_array
    q_mid = 0.5 * (q_full[:,:-1,:,:] + q_full[:,1:,:,:])
    u_mid = 0.5 * (u_full[:,:-1,:,:] + u_full[:,1:,:,:])
    v_mid = 0.5 * (v_full[:,:-1,:,:] + v_full[:,1:,:,:])

    # Calculate pressure jump
    dp = p_full[:,1:,:,:] - p_full[:,:-1,:,:]

    # Determine the fluxes and states
    fx = u_mid * q_mid * dp / g  # eastward atmospheric moisture flux (kg*m-1*s-1)
    fy = v_mid * q_mid * dp / g  # northward atmospheric moisture flux (#kg*m-1*s-1)
    cwv = q_mid * dp / g * a_gridcell / density_water  # column water vapor (m3)

    # convert fa_e, fa_n and cwv to xarray # TODO: move to preprocessing function
    def fs_to_xarray(fs, dimension_data, mid_levels, delta_pressure):
        fs_xr = xr.DataArray(
            data = fs,
            dims = ["time", "level", "latitude", "longitude"],
            coords = dict(
                time = dimension_data.time,
                latitude = dimension_data.latitude,
                longitude = dimension_data.longitude,
                mid_level = (["time", "level", "latitude", "longitude"], mid_levels),
                delta_pressure = (["time", "level", "latitude", "longitude"], delta_pressure)
            ) )
        return fs_xr

    fx = fs_to_xarray(fa_e, q, midpoints, dp)
    fy = fs_to_xarray(fa_n, q, midpoints, dp)
    cwv = fs_to_xarray(cwv, q, midpoints, dp)

    # Integrate fluxes and state
    upper_layer = above_boundary[:,:-1,:,:]
    lower_layer = ~upper_layer

    fx_lower = fx.where(lower_layer).sum(dim="level")  # kg*m-1*s-1
    fy_lower = fy.where(lower_layer).sum(dim="level")  # kg*m-1*s-1
    s_lower = cwv.where(lower_layer).sum(dim="level")  # m3

    fx_upper = fx.where(upper_layer).sum(dim="level")  # kg*m-1*s-1
    fy_upper = fy.where(upper_layer).sum(dim="level")  # kg*m-1*s-1
    s_upper = cwv.where(upper_layer).sum(dim="level")  # m3

    print(
        "Check calculation water vapor, this value should be zero:",
        (cwv.sum(dim="level") - (s_upper + s_lower)).sum().values,
    )

    tcwm3 = tcw * a_gridcell[np.newaxis, :] / density_water  # m3

    # Change units to m3, based on target frequency (not incoming frequency!)
    target_freq = config["target_frequency"]
    total_seconds = pd.Timedelta(target_freq).total_seconds()
    fx_upper *= total_seconds * (l_ew_gridcell / density_water)
    fx_lower *= total_seconds * (l_ew_gridcell / density_water)
    fy_upper *= total_seconds * (l_mid_gridcell[None, :, None] / density_water)
    fy_lower *= total_seconds * (l_mid_gridcell[None, :, None] / density_water)

    # Increased time resolution; states at midpoints, fluxes at the edges
    old_time = s_upper.time.values
    newtime_states = pd.date_range(old_time[0], old_time[-1], freq=target_freq)
    newtime_fluxes = newtime_states[:-1] + pd.Timedelta(target_freq) / 2

    s_upper = s_upper.interp(time=newtime_states).values
    s_lower = s_lower.interp(time=newtime_states).values
    fx_upper = fx_upper.interp(time=newtime_fluxes).values
    fy_upper = fy_upper.interp(time=newtime_fluxes).values
    fx_lower = fx_lower.interp(time=newtime_fluxes).values
    fy_lower = fy_lower.interp(time=newtime_fluxes).values
    precip = (precip.reindex(time=newtime_fluxes, method="bfill") / 4).values
    evap = (evap.reindex(time=newtime_fluxes, method="bfill") / 4).values

    # Stabilize horizontal fluxes
    fx_upper, fy_upper = get_stable_fluxes(fx_upper, fy_upper, s_upper)
    fx_lower, fy_lower = get_stable_fluxes(fx_lower, fy_lower, s_lower)

    # Determine the vertical moisture flux
    f_vert = get_vertical_transport(
        fx_upper,
        fx_lower,
        fy_upper,
        fy_lower,
        evap,
        precip,
        s_upper,
        s_lower,
        config["periodic_boundary"],
        config["kvf"],
    )

    # Save preprocessed data
    # Note: fluxes (dim: time) are at the edges of the timesteps,
    # while states (dim: time2) are at the midpoints and include next midnight
    # so the first state from day 2 will overlap with the last flux from day 1
    filename = f"{date.strftime('%Y-%m-%d')}_fluxes_storages.nc"
    output_path = output_dir / filename

    xr.Dataset(
        {  # TODO: would be nice to add coordinates and units as well
            "fx_upper": (
                ["time_fluxes", "lat", "lon"],
                fx_upper,
                {"units": "m**3"},
            ),
            "fy_upper": (
                ["time_fluxes", "lat", "lon"],
                fy_upper,
                {"units": "m**3"},
            ),
            "fx_lower": (
                ["time_fluxes", "lat", "lon"],
                fx_lower,
                {"units": "m**3"},
            ),
            "fy_lower": (
                ["time_fluxes", "lat", "lon"],
                fy_lower,
                {"units": "m**3"},
            ),
            "s_upper": (["time_states", "lat", "lon"], s_upper, {"units": "m**3"}),
            "s_lower": (["time_states", "lat", "lon"], s_lower, {"units": "m**3"}),
            "f_vert": (["time_fluxes", "lat", "lon"], f_vert, {"units": "m**3"}),
            "evap": (["time_fluxes", "lat", "lon"], evap, {"units": "m**3"}),
            "precip": (["time_fluxes", "lat", "lon"], precip, {"units": "m**3"}),
        },
        coords={
            "time_fluxes": newtime_fluxes,
            "time_states": newtime_states,
            "lat": u.latitude.values,
            "lon": u.longitude.values,
        },
    ).to_netcdf(output_path)
