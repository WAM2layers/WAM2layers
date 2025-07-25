
(opendap)=
# Using remote input data

WAM2layers now supports using remote input via OpenDAP. 40 (??) years of pre-processed ERA5 data are available at the 4TU data repository.

Below is an example configuration file that reads this input data:

```
# General
preprocessed_data_folder: https://opendap.4tu.nl/thredds/dodsC/data2/djht/00f7fa45-899e-4573-ae23-234f6c5193d0/1
output_folder: ./output_data
tracking_direction: backward

# Preprocessing (not used in this example but valid types still required)
filename_template: ./ERA5_{year}-{month:02d}-{day:02d}{levtype}_{variable}.nc
preprocess_start_date: "2021-07-01T00:00"
preprocess_end_date: "2021-07-15T00:00"
level_type: model_levels
levels: [20,40,60,80,90,95,100,105,110,115,120,123,125,128,130,131,132,133,134,135,136,137]

# Tracking
# tracking_domain: null  # use full domain of preprocessed data
tracking_domain: [0, 50, 10, 55]  # use subdomain of preprocessed data; west south east north
tracking_start_date: "1953-01-25T00:00"
tracking_end_date: "1953-02-02T00:00"

# Tagging
# tagging_region: tests/test_data/region_rhine.nc  # path must exist
tagging_region: [-4.0, 49.5, 6.0, 56.0]  # alternative to mask file; west south east north
tagging_start_date: "1953-01-31T00:00"
tagging_end_date: "1953-02-02T00:00"

input_frequency: "1h"
timestep: 600  # timestep in seconds
output_frequency: "1D"
restart: False
periodic_boundary: false
kvf: 3

```

Notice that the `preprocessed_data_folder` now points to the remote URL rather than a local path.

To use this configuration file, save it as `watersnoodramp.yaml` and start an experiment with

```
wam2layers track watersnoodramp.yaml
```

If you use this data for experiments, please cite the dataset as follows: ...
