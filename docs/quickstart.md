# How to use

## Quickstart
We recommend to create an editable installation in a dedicated conda environment (assuming you have an installation of Anaconda/Miniconda already):

```
# Create a dedicated conda environment
conda create -n wam2layers -c conda-forge python=3.9 jupyterlab cartopy matplotlib cmocean
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

A typical tracking experiment involves plotting and interactive
exploration. In fact, some of the plotting examples we provide are using
cartopy, which cannot be installed with pip. Thus, for an optimal experience we
recommend creating a dedicated conda environment for your moisture tracking
experiments, for example:

```
# Create a clean conda environment
conda create --name wamenv -c conda-forge python=3.9 jupyterlab cartopy matplotlib

# Activate the environment
conda activate wamenv
```

Now, you have two main options:

1. You do not want to modify the source code and use a tested version: install the latest 'release' of WAM2layers from PyPI (which should in principle coincide with the latest release on Zenodo):
```
# Install WAM2layers inside the new conda environment
pip install wam2layers
```

2. You want to modify the source code and create an editable installation from the latest master branch on Github (may be ahead of the PyPI/Zenodo release):
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

If you do not want to do anything interactive, the core of WAM2layers - the tracking functionality - is relatively
lightweight, and can simply be installed with pip:

```
# pip install git@github.com:WAM2layers/WAM2layers.git
or
pip install wam2layers
```

## Functionalities
WAM2layers consists of three basic functionalities: preprocessing, tracking and visualizing. All of them require a configuration file ('config_file'.yaml in the examples below). The visualize option just provides some basic examples of how you could visualize your data, but please modify that to your own needs. You can use these in the following way

preprocessing 
TODO: explain what happens (what goes in and what comes out, which files are used)
```
wam2layers preprocess 'config_filename'.yaml
```

You can visualize the output of the preprocessing, which is the input to the tracking by
```
wam2layers visualize input 'config_filename'.yaml
```

tracking (the core of WAM2layers):
TODO: explain what happens (what goes in and what comes out, which files are used)
```
wam2layers backtrack 'config_filename'.yaml
```

You can visualize the output of the preprocessing by
```
wam2layers visualize output 'config_filename'.yaml
```

TODO: explain the visualize both and snapshots options

## Download example data and configuration file

Before you can start tracking, you need some data. Input data for a simple example
case in which the pre-processing is already done is currently available through 
https://rebrand.ly/wam2layers-example-data. Once
stabilized, we will move it to a more persistent place, such as Zenodo.

The example data comes with a configuration file (.yaml file). It is best practice to save your configuration file outside of the WAM2layers repository.

You will need to update the paths in that file such that they point to your the location
where you stored the downloaded example data. TODO: 

For more information, see the documentation on [Case configuration](./config.md) and [Pre-processing](./preprocessing.md).

## Start an experiment

Now you are ready to start tracking! There are a few ways to use WAM2layers. 

One way is to use it from the command line:
```
wam2layers backtrack floodcase_202107.yaml
wam2layers visualize output floodcase_202107.yaml
```

This workflow is particularly useful if you are running experiments on HPC
systems.

Alternatively, if you prefer to work in an interactive Python environment such
as Jupyter Lab, you can import WAM2layers functionality and build your own scripts.

TODO: explain how it can by used with JupyterLab

```python
import wam2layers
config_file = "path_to/example_data/floodcase_2021.yaml"
wam2layers.tracking.backtrack.run_experiment(config_file)
```

## Running WAM2layers for your own case and data

### theory
Before using WAM2layers for your own case it is highly recommended to the [Theory](./theory.md) of tracking in WAM2layers 

### own case
Use the page [Case configuration](./config.md) to learn more about case configuration settings

### pre-processing
The default pre-processing works for ERA5 data for which example download scripts are provided in TODO:link to files. In case you just want to know more about what is happening during the pre-processing please refer to [Pre-processing](./prepocessing.md). 

### code adaptation and contribution
Possibly you want to adapt WAM2layers to track other climate data (e.g., output of RCMs and GCMs). In principle you can do so by only changing the Pre-processing and Case configuration without touching the core tracking module. You are encouraged to contribute your download scripts as well as prepocessing files to the WAM2layers Github repository.

## Code development
In case you are part of the development team (or want to become part of it), the [Developer's guide](./develop.md) is meant to help maintaining the code on Github, Zenodo and this documentation.



