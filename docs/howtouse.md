# How to use
TODO: introduce how this how to use works

## Installation instructions

A typical tracking experiment involves plotting and interactive
exploration. In fact, some of the plotting examples we provide are using
cartopy, which cannot be installed with pip. Thus, for an optimal experience we
recommend creating a dedicated conda environment (assuming you have an installation of Anaconda/Miniconda already) for your moisture tracking
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
pip install --editable .[develop]
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

## Running WAM2layers for your own case and data
TODO: merge the text below with the example?

### own case
Use the page [Case configuration](./config.md) to learn more about case configuration settings

### pre-processing
The default pre-processing works for ERA5 data for which example download scripts are provided in TODO:link to files. In case you just want to know more about what is happening during the pre-processing please refer to [Pre-processing](./prepocessing.md). 

### code adaptation and contribution
Possibly you want to adapt WAM2layers to track other climate data (e.g., output of RCMs and GCMs). In principle you can do so by only changing the Pre-processing and Case configuration without touching the core tracking module. You are encouraged to contribute your download scripts as well as prepocessing files to the WAM2layers Github repository.

## Code development
In case you are part of the development team (or want to become part of it), the [Developer's guide](./develop.md) is meant to help maintaining the code on Github, Zenodo and this documentation.