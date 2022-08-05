#!/bin/bash 

getSingleLevel=false
getModelLevel=true

# Single level data can be downloaded using era5cli
if [ "$getSingleLevel" = true ]
then
    era5cli hourly --startyear 2021 --months 7 --days 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 --area 60 -50 30 30 --variables total_precipitation evaporation surface_pressure 2m_temperature 2m_dewpoint_temperature 10m_u_component_of_wind 10m_v_component_of_wind 
fi

# Model level data have to be downloaded seperately
if [ "$getModelLevel" = true ]
then
    python get_modeldata.py
fi
