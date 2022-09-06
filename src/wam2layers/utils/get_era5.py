import cdsapi
import yaml
import pandas as pd
from pathlib import Path

c = cdsapi.Client()

# Read case configuration
from wam2layers.preprocessing.era5 import parse_config
config=parse_config("cases/era5_2021.yaml")

#Construct level list
if((config['levels'] == 'All') & (config['level_type'] == 'model_levels')):
    levlist = '1/to/137'
elif((config['levels'] == 'All') & (config['level_type'] == 'pressure_levels')):
    levlist = '1/2/3/5/7/10/20/30/50/70/100/125/150/175/200/225/250/300/350/400/450/500/550/600/650/700/750/775/800/825/850/875/900/925/950/975/1000'
else:
    s = '/'
    levlist = s.join(str(e) for e in config['levels'])

# Construct list of downloading times
if (config['data_freq'] == 'hourly'):
    timlist  = '00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00'
elif(config['data_freq'] == '3hourly'): 
    timlist = '00:00:00/03:00:00/06:00:00/09:00:00/12:00:00/15:00:00/18:00:00/21:00:00'
elif(config['data_freq'] == '6hourly'): 
    timlist = '00:00:00/06:00:00/12:00:00/18:00:00'
    
area = config['data_area'] # North, West, South, East. Default: global
grid = [config['data_resolution'], config['data_resolution']]  

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
        dtstr = date.strftime("%Y-%m-%d")

        # Create yearly parent dir if it doesn't exist yet
        parent = Path(config["input_folder"]).expand_user() / str(date.year)
        parent.mkdir(exist_ok=True, parents=True)
        
        # Create monthly subdirectories
        child = parent / str(date.month)
        child.mkdir(exist_ok=True)    # note: expanduser not necessary, already done above
        #Loop over model level variables
        for i in range(len(var_names)):
            outfile =  'ERA5_' + dtstr + "_" + var_names[i] + '_ml.nc'
            outfile = data_dir / outfile
            
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
             'grid': grid,       
        }, outfile) 
        
        #Loop over single level variables
        for i in range(len(sfc_long)):
            outfile = 'ERA5_' + dtstr + "_" + sfc_short[i] + '.nc' 
            outfile = data_dir / outfile
            
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
        dtstr = date.strftime("%Y-%m-%d")            
        
        #For first day of month, datadir should be created
        data_dir = Path(config["input_folder"] + "/" + str(date.year) ).expanduser()
        data_dir.mkdir(exist_ok=True, parents=True)
        
        #Subdirectory for every month 
        data_dir = Path(config["input_folder"] + "/" + str(date.year) + '/' + str(date.month) ).expanduser()
        data_dir.mkdir(exist_ok=True, parents=True)
        
        for i in range(len(varp)):
            outfile = 'ERA5_' + dtstr + "_" + varpshort[i] + '_pl_.nc'
            outfile = data_dir / outfile

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
            outfile =  'ERA5_' + dtstr + "_" + sfc_short[i] + '.nc' 
            outfile = data_dir / outfile

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