
(input-data)=
# Obtaining input data

Every tracking experiment requires data, which you might need to download.

**Meteorogical data**

Typically you will need wind speed in both directions and humidity at multiple
(pressure) levels and evaporation and precipitation. If the data comes at model
levels, you also need the corresponding pressure.

At this stage, the format of the data is not very strict yet. In the
[pre-processing step](./preprocess) we make sure that the data is ready for
tracking. However, if you want to use the pre-processing scripts that come
bundled with WAM2layers out-of-the-box, you need to make sure the data follows
the same structure as the [example data](example-data).

**Tracking region**

In addition to the time-dependent forcing data, WAM2layers also needs a
netcdf-file containing a variable called `tagging_region`, with values between 0
and 1. It must have the same lat/lon coordinates as the other input data. This
file is used to define the tagging region from where to track moisture.

(example-data)=
## Example data

To help you get started, we provide two example cases: 1) a backward tracking
case for an extreme precipitation event over the Eifel region in July 2021 and
2) a forward tracking case of tracking evaporation over the Volta region in
Ghana for July and August 1998. ERA5 data for these regions/periods are
available on 4TU and example configuration files are shipped with the package.

## File structure

A good file naming system makes it easy to process data for different date
ranges, and to automatically read data into the model. For example, for ERA5, we
are using the following filenaming structure:

```
path/to/data/<year>/<month>/ERA5_<year>_<month>_<day><levtype>_<variable>.nc

# for example:
path/to/data/2021/07/ERA5_2021_07_15_ml_q.nc
```

Here, the `<levtype>` parameter can be `_ml` for model levels, `_pl` for
pressure levels, and empty for surface variables. The `variable` name
corresponds to the name of the variable inside the netcdf file, so if our file
is called `......tp.nc` we expect to find a variable called `tp` inside this
file.

## Download scripts

Though it is not our main focus to provide data download utilities, we think it
may be helpful to share some commonly used scripts so that not everyone will
have to reinvent the wheel. Our download scripts for ERA5 (model and pressure
level) data are shipped with the repository and can be found in the [scripts
folder](https://github.com/WAM2layers/WAM2layers/tree/master/scripts). These
scripts are currently not considered core functionality of WAM2layers and are
therefore not part of the command line interface. Please see the docstrings
inside these scripts for instructions on how to use them.

```{note}
It would be nice to collect more data download scripts over time. If you want to
share your datascript, we will be happy to help you in adding it to our
collection. Please see the [contributing](../develop.md) chapter of this
documentation.
```
