import cdsapi
import yaml
import pandas as pd
c = cdsapi.Client()

# Read case configuration
with open("cases/era5_2021.yaml") as f:
    config = yaml.safe_load(f)
    
datelist = pd.date_range(
    start=config["preprocess_start_date"], end=config["preprocess_end_date"], freq="d"
)

if(config['levels'] == 'All'):
    levlist = '1/to/137'
else:
    s = '/'
    levlist = s.join(str(e) for e in config['levels'])
    
timlist  = '00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00'
area = [60, -50,30,30] # North, West, South, East.          Default: global
grid = [0.25, 0.25]  # Latitude/longitude grid.           Default: 0.25 x 0.25

######################################################################
# For model levels the following variables are downloaded:
# u,v, q on selected levels
# tp (total precipitation), e (evaporation) and sp (surface pressure)
######################################################################
if(config['level_type'] == 'model_levels'):
    var_names = ["u","v","q"]
    var_ids = [131,132,133]
    
    sfc_short = ['tp','e','sp']
    sfc_long = ['mean_sea_level_pressure', 'evaporation', 'total_precipitation']
    
    for date in datelist:
        #For first day of month, datadir should be created
        dtstr = date.strftime("%Y-%m-%d")
        
        #Loop over model level variables
        for i in range(len(var_names)):
            outfile = config['input_folder'] + '/ERA5_' + dtstr + "_" + var_names[i] + '_ml_.nc'
            
            c.retrieve('reanalysis-era5-complete', {
            'class': 'ea',
            'date': dtstr,
            'expver': '1',
            'levelist': levlist,
            'levtype': 'ml',
            'param': var_ids[i],  
            'stream': 'oper',
            'time': timlist,
            'type': 'an',
            'format': 'netcdf',
             'area': area,
             'grid': grid,       # Latitude/longitude grid.           Default: 0.25 x 0.25
        }, outfile) 
        
        #Loop over single level variables
        for i in range(len(sfc_long)):
            outfile = config['input_folder'] + '/ERA5_' + dtstr + "_" + sfc_short[i] + '.nc'  
            c.retrieve(
            'reanalysis-era5-single-levels',
            {
                'product_type': 'reanalysis',
                'variable': sfc_long[i],
                'date': dtstr,
                'time': timlist,
                'area': area, 
                'grid': grid,      
                'format': 'netcdf',
            },
            outfile)
######################################################################
# For pressure levels the following variables are downloaded:
# u,v, q on selected pressure levels
# tp (total precipitation), e (evaporation) and sp (surface pressure)
# + t2m (2 meter temperature), 
######################################################################
elif(config['level_type'] == 'pressure_levels'):
    print("Getting pressure levels")
    
    varp= ['specific_humidity', 'u_component_of_wind', 'v_component_of_wind']  
    varpshort = ['q','u','v']
    
    sfc_short = ['tp','e','sp','d2m','u10','v10']
    sfc_long = ['mean_sea_level_pressure', 'evaporation', 'total_precipitation','2m_dewpoint_temperature','10m_u_component_of_wind', '10m_v_component_of_wind']
     
    for date in datelist:
        #For first day of month, datadir should be created
        dtstr = date.strftime("%Y-%m-%d")            
        for i in range(len(varp)):
            outfile = config['input_folder'] + '/ERA5_' + dtstr + "_" + varpshort[i] + '_pl_.nc'
            c.retrieve(
            'reanalysis-era5-pressure-levels',
            {
                'product_type': 'reanalysis',
                'variable': varp[i],
                'pressure_level': levlist,
                'date': dtstr,
                'time': timlist,
                'area': area, 
                'grid': grid, 
                'format': 'netcdf',
            },
            outfile)
            
            #Loop over single level variables
        for i in range(len(sfc_long)):
            outfile = config['input_folder'] + '/ERA5_' + dtstr + "_" + sfc_short[i] + '.nc'  
            c.retrieve(
            'reanalysis-era5-single-levels',
            {
                'product_type': 'reanalysis',
                'variable': sfc_long[i],
                'date': dtstr,
                'time': timlist,
                'area': area,
                'grid': grid,  
                'format': 'netcdf',
            },
            outfile)