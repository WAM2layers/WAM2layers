# User guide

In this section we will guide you through an example of how we are using
WAM2layers in practice. By the end of this page, you will be familiar with the
different functionalities of WAM2layers, and have a basic understanding of how
the program works. Moreover, you will understand how you can set up your own
experiments with WAM2layers. This guide can also be used as a reference to
quickly find some specific information.

A typical experiment may consist of the following steps:

1. [Obtain input data](input-data)
1. [Pre-process the data](pre-processing)
1. [Perform the actual tracking](tracking)
1. [Analyse the outputs](analysis)

Not all of these steps require WAM2layers. We describe them nonetheless as they
provide a lot of relevant context on how WAM2layers is used in practice.


(cli)=
## WAM2layers' command line interface

WAM2layers is written in Python, but unlike many popular Python packages, it is
not primarily intended as a library that you can import like numpy. While it is
possible to import and [use wam2layers in an interactive session](interactive), this is
probably not the best experience at the moment. Instead, we offer a command line
interace (CLI), much like pip, conda, git, cdo, et cetera. It looks like this:

```
# Get help
wam2layers --help
wam2layers track --help

# Display the version of wam2layers you are using:
wam2layers --version

# Preprocess data (you need to have downloaded the data)
wam2layers preprocess era5 floodcase_202107.yaml

# Perform tracking
wam2layers track floodcase_202107.yaml

# Make some basic plots
wam2layers visualize output floodcase_202107.yaml
```

All these commands are explained in detail below.

```{admonition} From command line to python code

Usually you will not need to bother with this, but it can be nice to know how
the command line statements are connected to the source code. Three things are key
here:

1. In essence, WAM2layers is just a collection of scripts.
1. The file
   [`pyproject.toml`](https://github.com/WAM2layers/WAM2layers/blob/main/pyproject.toml),
   which is used during installation, tells your system that the command
   `wam2layers` should point to the file `src/wam2layers/cli.py`.
1. We use [`click`](https://click.palletsprojects.com/en/8.1.x/) to further
   process everything that comes after the `wam2layers` command. For example,
   `wam2layers track` is redirected to the appropriate module
   (`src/wam2layers/tracking/backtrack.py:run_experiment` or
   `src/wam2layers/tracking/forwardtrack.py:run_experiment`). This function takes the
   configuration file as input and passes it on to a function called
   `run_experiment`.
```

The command line interface is very useful when you want to do multiple
experiments with different settings. You only ever edit the config file, not the
source code. It is also useful for running the model on an HPC system. Moreover,
the command `wam2layers` can be used from anywhere on your system. However, if
you want, you can also run

```
# Call the backtrack script directly with python. Requires the full path.
python src/wam2layers/tracking/backtrack.py config-file.yaml
```

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
collection. Please see the [contributing](./develop.md) chapter of this
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

(pre-processing)=
## Pre-processing

During the pre-processing step we make sure to convert the raw data into
something that the tracking script understands. This typically involves
accumulating the data into two layers, deriving moisture fluxes (which are
typically not included in standard model output) from variables that are
commonly available, and converting the data to the right units and a regular
grid if that isn't the case already.

### Built-in preprocessing functionality

WAM2layers comes with built-in preprocessing functionality for ERA5 data. You can
use the following command:

```
# Preprocess era5 data using the wam2layers builting functionality
wam2layers preprocess era5 config_file.yaml
```

where `config_file.yaml` is the path to your configuration file. This file
should have settings on the date range for which you want to run the
preprocessing, and also about the location where the raw data are stored and
what filename pattern they follow. For more information, see [](./config) or
have a look at the example config file
[here](https://github.com/WAM2layers/WAM2layers/blob/main/example-config.yaml).

```{note}
For now WAM2layers only contains preprocessing code for era5. We think it would
be nice to add preprocessing functionality for more datasets over time.
```

### Preprocessing other datasets

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
very grateful if you are willing to [contribute your code](./develop.md) so
others can benefit from it as well.
```

### Data checking utilities

To aid in pre-processing, WAM2layers comes with some basic data checking
utilities. Specifically, the function `check_input` takes a single input file
(opened with `xarray`) as input and checks it against most of the requirements
enumerated above. You can use this to get some reassurance (or detect issues).


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

(tracking)=
## Tracking

The core of WAM2layers are the tracking routines. The model includes forward and
backward tracking. Time tracking, distance tracking, and moisture recycling can
be added in future updates.

Assuming you have a preprocessed dataset and prepared a configuration file for
your experiment, running WAM2layers is as simple as:

```
wam2layers track config-file.yaml
```

where `config-file.yaml` is the path to your configuration file. Among others,
this file should have settings on the tracking direction (forward/backward),
the date range for which you want to run track, and also about the location where the
preprocessed data are stored and where the output will be stored.
For more information, see [](./config) or have a look at the example config file
[here](https://github.com/WAM2layers/WAM2layers/blob/main/example_config.yaml).

```{tip}
The actual code that is executed can be found in
`src/wam2layers/tracking/backtrack.py` or `src/wam2layers/tracking/forwardtrack.py`. This script reads the configuration,
loads the preprocessed data step by step, calculates the proportion of tracked
moisture in each grid cell, and writes the output to a file. You can configure
how often output is written.

For more information on the equations being solved, see the
[theory](./theory.md) chapter. If you want to know how the command line call
ends up in this python script, see the section on [WAM2layers' command line
interface](cli)
```

(analysis)=
## Analysis and visualization

It can be useful to quickly make some standard plots of the input and output.
To this end, WAM2layers includes a few visualization scripts.

```
# See what's available
wam2layers visualize --help

# Visualize the pre-processed data
wam2layers visualize input config-file.yaml

# Visualize the output data after tracking
wam2layers visualize output config-file.yaml

TODO: explain the visualize both and snapshots options
```

The configuration file is used to find the path to the data.

Of course, you are free to analyse and plot the data completely to your own
liking. These utilities are just meant to quickly inspect the data. If you find
yourself constantly making another type of plot, and you think it might be
useful to have this as a standard plot in WAM2layers, please [reach out to us
through GitHub](https://github.com/waM2layers/waM2layers/issues/new).

(interactive)=
## Using WAM2layers in an interactive session

While WAM2layers was designed to be run from the command line, it is possible to
run it in an interactive session such as IPython or Jupyter Lab. You can import
WAM2layers functionality and build your own scripts. For example:

```python
import wam2layers
config_file = "path_to/example_data/floodcase_2021.yaml"
wam2layers.tracking.backtrack.run_experiment(config_file)
```

You can also import specific functions. However, you should be conscious that
the code was not designed to be used in this way, so it might not be the most
intuitive, and you might run into unexpected behaviour. Moreover, in maintaining
the model, we cannot guarentee backward compatibility for this type of use. So
use it at your own risk.

Perhaps, at some point, it would be nice to build a more user-friendly Python
API, for example following the ["basic model
interface"](https://bmi.readthedocs.io/en/stable/):

```python
from wam2layers import BacktrackModel

config_file = "path_to/example_data/floodcase_2021.yaml"

model = BacktrackModel()
model.initialize(config_file)

while model.time < model.end_time:
    model.update()

model.get_value('s_upper')
```
