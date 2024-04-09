
(input-data)=
## Input data

Every tracking experiment requires data, which you might need to download.
Typically you will need wind speed in both directions and humidity at multiple
(pressure) levels and evaporation and precipitation. If the data comes at model
levels, you also need the corresponding pressure.

### Download scripts

Though it is not our main focus to provide data download utilities, we think it
may be helpful to share some commonly used scripts so that not everyone will
have to reimplement the wheel. Our download scripts for ERA5 (model and pressure
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

### File structure

A good file naming system (a.k.a. *data reference syntax*) makes it easy to
process data for different date ranges, and to automatically read data inside
our preprocessing scripts. In our download scripts we are using the following
filenaming structure:

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


### Tracking region

In addition to the time-dependent forcing data, WAM2layers also needs a netcdf-file containing a variable called `tagging_region`, with values between 0 and 1. It must have the same lat/lon coordinates as the other
input data. This file is used to define the tagging region from where to track
moisture.

### Example data

Input data for a simple example
case in which the pre-processing is already done is currently available through
https://rebrand.ly/wam2layers-example-data.

In case you're struggling with the steps above and want to have a quick
impression of how WAM2layers works, we have made a small, pre-processed data
available at <https://rebrand.ly/wam2layers-example-data>. Once stabilized, we
will move it to a more persistent place, such as Zenodo. You can use this
dataset and proceed with the next steps.
