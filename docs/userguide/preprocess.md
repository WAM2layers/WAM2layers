
(pre-processing)=
# Pre-processing

During the pre-processing step we make sure to convert the raw data into
something that the tracking script understands. This typically involves
accumulating the data into two layers, deriving moisture fluxes (which are
typically not included in standard model output) from variables that are
commonly available, and converting the data to the right units and a regular
grid if that isn't the case already.

## Built-in preprocessing functionality

WAM2layers comes with built-in preprocessing functionality for ERA5, ARCO-ERA5 and CMIP
data. You can use the following commands:

```
# Preprocess data using the wam2layers built-in functionality
wam2layers preprocess era5 config_file.yaml
wam2layers preprocess arco-era5 config_file.yaml
wam2layers preprocess cmip config_file.yaml
```

where `config_file.yaml` is the path to your configuration file. This file
should have settings on the date range for which you want to run the
preprocessing, and also about the location where the raw data are stored and
what filename pattern they follow. For more information, see [](./config) or
have a look at the example config file
[here](https://github.com/WAM2layers/WAM2layers/blob/main/example-config.yaml).

### ARCO-ERA5

[ARCO-ERA5](https://github.com/google-research/arco-era5) is a version of ERA5 hosted
publicly by Google.
When you use this input data for the preprocessor, you don't need to download any files
to your computer. The input data will be pulled in on-the-fly.

However, there are some limitation to this. Due to the way the dataset is stored, all
latitude/longitude values, as well as all pressure levels are retrieved at a time.
You cannot download subsets of this data.

For this reason the preprocessor will *always download data for the entire globe*.
The pressure levels you specify in your config will be kept, any others are still
downloaded but "thrown away". For this reason you could pre-process all pressure levels
but only use some in your analysis (if disk space allows).

```{note}
If you and other people in your group/institute make use of this data, it could be 
useful to store it somewhere where all of you can have shared access, and preprocess
the entire globe/all pressure levels once.
```

## Preprocessing other datasets

If you want to use another dataset, you need to make sure that it follows the
same standards. To give you an impression, here is an example structure of a
preprocessed netCDF file:

```
Dimensions:    (time: 25, longitude: 321, latitude: 121)
Coordinates:
  * time       (time) datetime64[ns] 2021-07-01 ... 2021-07-02
  * longitude  (longitude) float32 -50.0 -49.75 -49.5 -49.25 ... 29.5 29.75 30.0
  * latitude   (latitude) float32 60.0 59.75 59.5 59.25 ... 30.5 30.25 30.0
Data variables:
    fx_upper   (time, latitude, longitude) float32 ...
    fy_upper   (time, latitude, longitude) float32 ...
    fx_lower   (time, latitude, longitude) float32 ...
    fy_lower   (time, latitude, longitude) float32 ...
    s_upper    (time, latitude, longitude) float32 ...
    s_lower    (time, latitude, longitude) float32 ...
    evap       (time, latitude, longitude) float32 ...
    precip     (time, latitude, longitude) float32 ...
```

This pre-processed dataset adheres to the following requirements:

- Data should be stored in netcdf files, one file per day.
- Precipitation and evaporation should both be positive
- The states (`s_upper` and `s_lower`) should be positive as well
- States are given in units of "kg m-2", fluxes in "kg m-1 s-1"
- `evap` and `precip` are given in units of "kg m-2 s-2"
- Latitude should be decreasing, time and longitude increasing.

```{note}
If you need help in pre-processing your data, please don't hesitate to reach
out, for example through
[GitHub](https://github.com/WAM2layers/WAM2layers/issues/new). We would also be
very grateful if you are willing to [contribute your code](../develop.md) so
others can benefit from it as well.
```


<!-- TODO: update this and make it work and look nice(r) -->
<!-- ## Data checking utilities

To aid in pre-processing, WAM2layers comes with some basic data checking
utilities. Specifically, the function `check_input` takes a single input file
(opened with `xarray`) as input and checks it against most of the requirements
enumerated above. You can use this to get some reassurance (or detect issues). -->
