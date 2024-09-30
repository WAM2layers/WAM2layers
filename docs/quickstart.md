# Quickstart

To help you get started, we provide two example cases:

- a backward tracking case for an extreme precipitation event over the Eifel
  region in July 2021
- a forward tracking case of tracking evaporation over the Volta region in Ghana
  for July and August 1998.

ERA5 data for these regions/periods are
available on 4TU and example configuration files are shipped with the package.

This quickstart gives a good impression of all the steps needed to run these
example cases.

```sh
# Optional: Create a new conda environment for WAM2layers and activate it
conda create --name wamenv python=3.11
conda activate wamenv

# Install wam2layers in your current python environment
pip install wam2layers

# Check the version of wam2layers
wam2layers --version

# See the available commands/options
wam2layers --help
wam2layers track --help

# Download example input data (pick one)
wam2layers download example-input-eiffel
wam2layers download example-input-volta

cd example-input-eiffel  # or example-input-volta

# Prepare the data for a tracking experiment
wam2layers preprocess era5 config-eiffel.yaml  # or example-config-volta

# Run the tracking experiment
wam2layers track config-eiffel.yaml  # or example-config-volta

# Make a default plot of the results
wam2layers visualize output config-eiffel.yaml  # or example-config-volta
```

A detailed explanation for each of these steps is available in the [user
guide](./userguide/index).
