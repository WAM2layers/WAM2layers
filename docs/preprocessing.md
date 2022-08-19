# Pre-processing

The tracking code makes several assumptions about the incoming data. It is the
responsibility of the user to make sure that the input data is correctly
pre-processed. Some example scripts for preprocessing ERA5 data are [shipped
with the
repository](https://github.com/WAM2layers/WAM2layers/tree/master/wam2layers/preprocessing).

The final input data into the tracking should have the following
characteristics:

- Data should be stored in netcdf files, one file per day.
- The file should contain the following variables:
   - fx_upper: eastward moisture transport in the upper layer
   - fy_upper: northward moisture transport in the upper layer
   - fx_lower: eastward moisture transport in the lower layer
   - fy_lower: northward moisture transport in the lower layer
   - w_upper: total water in the upper layer
   - w_lower: total water in the lower layer
   - evap: evaporation from the surface
   - precip: precipiation
- All variables should be in units of m3
- Precipitation and evaporation should both be positive.

Here is an example of a preprocessed netCDF file. Note that the latitude,
longitude, and time may vary for your data.

```
Dimensions:   (time: 95, lat: 121, lon: 321, time2: 96)
Dimensions without coordinates: time, lat, lon, time2
Data variables:
    fx_upper  (time, lat, lon) float64 ...
    fy_upper  (time, lat, lon) float64 ...
    fx_lower  (time, lat, lon) float64 ...
    fy_lower  (time, lat, lon) float64 ...
    w_upper   (time, lat, lon) float64 ...
    w_lower   (time, lat, lon) float64 ...
    evap      (time, lat, lon) float32 ...
    precip    (time, lat, lon) float32 ...
```
