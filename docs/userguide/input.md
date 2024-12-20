
(input-data)=
# Obtaining input data

Every tracking experiment requires data, which you might need to download.

* **Meteorogical data**

    Typically you will need wind speed in both directions and humidity at multiple
    (pressure) levels and evaporation and precipitation. If the data comes at model
    levels, you also need the corresponding pressure.


* **Tagging region**
```{eval-rst}
    In addition to the time-dependent forcing data, WAM2layers also needs a
    tagging region. This is used to define the tagging region from where to
    track moisture.
    The region can be described by a netCDF file, shapefile, or simple
    bounding box. For more info see :class:`wam2layers.config.Config.tagging_region`

    Alternatively, on the respository `some notebooks are available <https://github.com/WAM2layers/WAM2layers/tree/main/notebooks>`_ which tackle more advances tagging region cases.
```


At this stage, the format of the data is not very strict yet. In the
[pre-processing step](./preprocess) we make sure that the data is ready for
tracking. However, if you want to use the pre-processing scripts that come
bundled with WAM2layers out-of-the-box, you need to make sure the data follows
the same structure as the [example data](example-data).

(example-data)=
## Example data

To help you get started, we provide two example cases:

1. A backward tracking case for an extreme precipitation event over the Eifel
   region in July 2021.
2. A forward tracking case of evaporation over the Volta region in
   Ghana for July and August 1998.

ERA5 data for these regions/periods are available on 4TU
([Eiffel](https://doi.org/10.4121/f9572240-f179-4338-9e1b-82c5598529e2);
[Volta](https://doi.org/10.4121/bbe10a2a-39dc-4098-a69f-0f677d06ecdd))
and can easily be downloaded with WAM2layers:

```sh
wam2layers download example-input-eiffel
wam2layers download example-input-volta
```

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
