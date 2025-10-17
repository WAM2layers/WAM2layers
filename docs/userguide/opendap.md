(opendap)=

# Using Remote Preprocessed Data

WAM2layers supports loading preprocessed data directly from **OpenDAP**, allowing you to skip downloading and preprocessing steps and move straight to running tracking experiments.

As a basis, [Bakels et al. (2025)](https://doi.org/10.4121/00f7fa45-899e-4573-ae23-234f6c5193d0.v1) provides 83 years of preprocessed ERA5 data available via the 4TU ResearchData repository. This dataset:

* covers the period **1941–2024**
* is regridded to a resolution of **0.5° × 0.5°**
* uses the 'standard' model level configurations with 22 model levels

## Example Configuration

In the configuration file the preprocessed data folder uses a remote dataset:

```yaml
# General settings
preprocessed_data_folder: https://opendap.4tu.nl/thredds/dodsC/data2/djht/00f7fa45-899e-4573-ae23-234f6c5193d0/1
```

👉 **Note:** `preprocessed_data_folder` points to a **remote URL** instead of a local directory.


## Running the Experiment

Example configuration files that work with OpenDAP can be found [here](https://github.com/WAM2layers/WAM2layers/tree/main/configs). 
You run them as any other configuration file, for example:

```bash
wam2layers track opendap-config-volta.yaml
```


```{Admonition} Convenience vs performance
:class: tip

Fetching data on the fly is slower. In our test case, running with local files was about 5× faster (excluding the initial download time). Performance will vary with your internet connection, so weigh the convenience of remote access against speed when choosing your setup.
```

## Citation

If you use this dataset, please cite it as:

**Bakels, Lucie, van der Ent, R. J. (Ruud), Wang-Erlandsson, Lan, & Kalverla, Peter. (2025). *WAM2layers preprocessed data from 01-01-1941 until 31-10-2014*. 4TU.ResearchData. [https://doi.org/10.4121/00f7fa45-899e-4573-ae23-234f6c5193d0.v1](https://doi.org/10.4121/00f7fa45-899e-4573-ae23-234f6c5193d0.v1) (CC BY 4.0)**
