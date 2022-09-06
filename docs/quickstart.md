# How to use

## Quickstart
We recommend to create an editable installation in a dedicated conda environment:

```
# Create a dedicated conda environment
conda create -n wam2layers -c conda-forge python=3.9 jupyterlab cartopy matplotlib
conda activate wam2layers

# Clone the source code repository
git clone git@github.com:WAM2layers/WAM2layers.git
cd WAM2layers

# Install the package from source
pip install -e .

# Download example data and update the paths in the configuration file
https://rebrand.ly/wam2layers-example-data

# Run example case
wam2layers backtrack floodcase_202107.yaml
```

## Detailed installation instructions

The core of WAM2layers - the tracking functionality - is relatively
lightweight, and can simply be installed with pip:

```
pip install git@github.com:WAM2layers/WAM2layers.git
# pip install wam2layers  # will come soon
```

However, a typical tracking experiment also involves plotting and interactive
exploration. In fact, some of the plotting examples we provide are using
cartopy, which cannot be installed with pip. Thus, for an optimal experience we
recommend creating a dedicated conda environment for your moisture tracking
experiments, for example:

```
# Create a clean conda environment
conda create -n wam2layers -c conda-forge python=3.9 jupyterlab cartopy matplotlib

# Activate the environment
conda activate wam2layers

# Install WAM2layers inside the new conda environment
pip install wam2layers
```

If you want to modify the source code, you can create an editable installation:
```
# Clone the source code repository
git clone git@github.com:WAM2layers/WAM2layers.git

# Enter the source code directory
cd WAM2layers/

# Create an editable installation - with development dependencies
pip install -e .[develop]
```

Now, if you edit the source code of WAM2layers, the model will run with your
updated code.

## Download example data

Before we can start tracking, we need some data. Input data for a simple example
case is currently available via https://rebrand.ly/wam2layers-example-data. Once
stabilized, we will move it to a more persistent place, such as Zenodo.

You will see that the example data comes with a configuration file. You will
need to update the paths in that file such that they point to your the location
where you stored the downloaded example data.

For more information, see the documentation on configuration and pre-processing.

## Start an experiment

Now you are ready to start tracking! There are a few ways to use WAM2layers. One
way is to use it from the command line:
```
wam2layers backtrack floodcase_202107.yaml
wam2layers visualize output floodcase_202107.yaml
```

This workflow is particularly useful if you are running experiments on HPC
systems.

Alternatively, if you prefer to work in an interactive Python environment such
as Jupyter Lab, you can import WAM2layers functionality and build your own scripts.

```python
import wam2layers
config_file = "path_to/example_data/floodcase_2021.yaml"
wam2layers.tracking.backtrack.run_experiment(config_file)
```
