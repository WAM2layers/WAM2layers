# Pre-processing

The tracking code makes several assumptions about the incoming data. It is the
responsibility of the user to make sure that the input data is correctly
pre-processed. Some example scripts for preprocessing ERA5 data are [shipped
with the
repository](https://github.com/WAM2layers/WAM2layers/tree/master/wam2layers/preprocessing).

Here is an example of a preprocessed netCDF file. Note that the latitude,
longitude, and time may vary for your data.

```
Dimensions:    (time: 25, longitude: 321, latitude: 121)
Coordinates:
  * time       (time) datetime64[ns] 2021-07-01 ... 2021-07-02
  * longitude  (longitude) float32 -50.0 -49.75 -49.5 -49.25 ... 29.5 29.75 30.0
  * latitude   (latitude) float32 60.0 59.75 59.5 59.25 ... 30.5 30.25 30.0
Data variables:
    fx_upper   (time, latitude, longitude) float64 ...
    fy_upper   (time, latitude, longitude) float64 ...
    fx_lower   (time, latitude, longitude) float64 ...
    fy_lower   (time, latitude, longitude) float64 ...
    s_upper    (time, latitude, longitude) float64 ...
    s_lower    (time, latitude, longitude) float64 ...
    evap       (time, latitude, longitude) float32 ...
    precip     (time, latitude, longitude) float32 ...
```

This dataset adheres to the following requirements:

- Data should be stored in netcdf files, one file per day.
- Importantly, the midnight of the next day should also be included (except
  perhaps for the last day)
- Precipitation and evaporation should both be positive
- The states (`s_upper` and `s_lower`) should be positive as well
- States are given in units of "kg m-2", fluxes in "kg m-1 s-1"
- `evap` and `precip` are given in units of "kg m-2 s-2"
- Latitude should be decreasing, time and longitude increasing.

## Tracking region

In addition to the time-dependent forcing data, WAM2layers also needs a file
called `source_region.nc` containing a variable called `source_region`, with
values between 0 and 1. It must have the same lat/lon coordinates as the other
input data. This file is used to define the source region from where to track
moisture.

## Data checking utility functions

To aid in pre-processing, WAM2layers comes with some handy data checking
utilities. Specifically, the function `check_input` takes a single input file
(opened with `xarray`) as input and checks it against most of the requirements
enumerated above.
