(opendap)=

# Using Remote Input Data

WAM2layers supports loading input data directly from **OpenDAP**, allowing you to skip downloading and preprocessing steps and move straight to running tracking experiments.

As a basis, [Bakels et al. (2025)](https://doi.org/10.4121/00f7fa45-899e-4573-ae23-234f6c5193d0.v1) provides 73 years of preprocessed ERA5 data available via the 4TU ResearchData repository. This dataset:

* covers the period **1941–2014**
* is regridded to a resolution of **0.5° × 0.5°**

## Example Configuration

The example below shows how to set up a configuration file that uses the remote dataset:

```yaml
# General settings
preprocessed_data_folder: https://opendap.4tu.nl/thredds/dodsC/data2/djht/00f7fa45-899e-4573-ae23-234f6c5193d0/1
output_folder: ./output_data
tracking_direction: backward

# Preprocessing (not used here but required)
filename_template: ./ERA5_{year}-{month:02d}-{day:02d}{levtype}_{variable}.nc
preprocess_start_date: "2021-07-01T00:00"
preprocess_end_date: "2021-07-15T00:00"
level_type: model_levels
levels: [20,40,60,80,90,95,100,105,110,115,120,123,125,128,130,131,132,133,134,135,136,137]

# Tracking
tracking_domain: [0, 50, 10, 55]  # subdomain (west, south, east, north)
tracking_start_date: "1953-01-25T00:00"
tracking_end_date: "1953-02-02T00:00"

# Tagging
tagging_region: [-4.0, 49.5, 6.0, 56.0]
tagging_start_date: "1953-01-31T00:00"
tagging_end_date: "1953-02-02T00:00"

# Simulation settings
input_frequency: "1h"
timestep: 600
output_frequency: "1D"
restart: false
periodic_boundary: false
kvf: 3
```

👉 **Note:** `preprocessed_data_folder` points to a **remote URL** instead of a local directory.


## Running the Experiment

Save the configuration as `watersnoodramp.yaml` and run:

```bash
wam2layers track watersnoodramp.yaml
```


```{Admonition} Convenience vs performance
:class: tip

Fetching data on the fly is slower. In our test case, running with local files was about 5× faster (excluding the initial download time). Performance will vary with your internet connection, so weigh the convenience of remote access against speed when choosing your setup.
```

## Citation

If you use this dataset, please cite it as:

**Bakels, Lucie, van der Ent, R. J. (Ruud), Wang-Erlandsson, Lan, & Kalverla, Peter. (2025). *WAM2layers preprocessed data from 01-01-1941 until 31-10-2014*. 4TU.ResearchData. [https://doi.org/10.4121/00f7fa45-899e-4573-ae23-234f6c5193d0.v1](https://doi.org/10.4121/00f7fa45-899e-4573-ae23-234f6c5193d0.v1) (CC BY 4.0)**
