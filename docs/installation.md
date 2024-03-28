# Setup

This page describes how to install WAM2Layers for various use cases, from using
the model as is, to more advanced use cases for when you want to make changes to
the code. By the end of this walkthrough you should have a working installation
of WAM2Layers and understand how you ended up there.

If you want to make changes to the code and contribute them so that others can
also benefit from you work (highly apreciated!), you may want to follow the
instructions in the [contributing](./contributing) section of the documentation.

## Quickstart

If you know your way around, the following quick reference might be all you need:

Install WAM2layers in an existing Python environment:
```
pip install wam2layers[viz]
```

Create a new Jupyter Lab environment with conda and install wam2layers inside it.
```
conda create --name wamenv -c conda-forge python=3.11 jupyterlab cartopy matplotlib
conda activate wamenv
pip install wam2layers
```

If you're not completely sure what these commands do and how they work, please
read on. We'll explain everything in detail.

## Requirements

To run WAM2Layers you need a computer and some kind of command line, such as the
default terminal included in Mac or Linux. On Windows we recommend using Ubuntu
through [WSL](). If you're using existing infrastructure such as an managed
JupyterLab environment or a supercomputer you will also typically have access to
a terminal.

## Installation

The core of WAM2layers - the tracking functionality - is relatively lightweight.
It requires Python 3 and can be installed with pip (the default installer for
Python packages; often it comes included with Python). To check if Python and
pip are available on your system, you can use `which`:

```
# Check if python and pip are installed
which python3
which pip
```

If this returns nothing, you don't have Python/pip. In that case, we recommend 
installing Python with conda (see below).

If you do have Python and pip installed, you can install wam2layers with one of the 
following commands, depending on your requirements:

```
pip install wam2layers
pip install wam2layers[viz]  # includes packages required for visualization
pip install wam2layers[complete]  # include all optional packages (for development)
```

### Conda installation

If you are less experienced with setting up Python virtual environments,
it can be more convenient to use Conda.
Additionally, some packages might not be available from pip on all systems.

Again, you can use `which conda` to see if conda is available on your system. 
If it's not, [here](https://docs.anaconda.com/free/miniconda/index.html) is a 
really nice walkthrough for installing miniconda.

If you have conda installed, to set up a dedicated environment for your
moisture tracking experiments do:

```
# Create a fresh conda environment
conda create --name wamenv -c conda-forge python=3.11 jupyterlab cartopy matplotlib

# Activate the environment
conda activate wamenv
```

Now, if you simply want to use WAM2Layers as is, without modifying the source code,
you can install it from PyPI (the official python package repository):

```
# Install the latest 'release' from PyPI
pip install wam2layers

# You can also install a specific version:
pip install wam2layers=3.0.0-beta5
```

### Requirements for modifying the source code: git
In principle, the previous steps are enough to start using WAM2Layers.
In our experience, however, it is often useful and insightful to obtain a copy
of the source code and create an editable installation. This will allow you to
look at the source code to see what is going on inside the model. You can then
also make (experimental) changes to the model and see how it affects the
results.
The source code for WAM2Layers is hosted on 
GitHub. Creating a local copy requires git (tip: `which git`).

Beware that the code on GitHub might include recent updates that have
not yet been included in the latest official release. This can be an advantage,
but also a downside as it may be a bit more unstable.

If you want to contribute your changes, raise an issue, request new features, et
cetera, you may also want to sign up for a free GitHub account if you don't have
one already. For more information on contributing code, see the [contribution
guidelines](./contributing) chapter of the documentation.

First, obtain a copy of the source code from GitHub and enter the new directory:

```
# Clone the source code repository
git clone https://github.com/WAM2layers/WAM2layers.git

# Enter the source code directory
cd WAM2layers/
```

Local code can also be installed with pip, but instead of pointing to PyPI as
above, point to the current working directory (denoted by `.`). Additionally,
add the `--editable` flag, and `[complete]` to install all optional packages.

```
# Make an editable installation of the source code in this directory (.)
pip install --editable .[complete]
```

Now, if you edit the source code of WAM2layers, the model will run with your
updated code. If you want to learn more about how pip recognized that your
current working directory contained an installable package, have a look at [the
official
documentation](https://packaging.python.org/en/latest/tutorials/packaging-projects/).

## Next steps

That should be all. To confirm that WAM2layers is installed correctly, you can use `which`, or ask for help:

```
# See if WAM2layers is installed on your system
which wam2layers

# See which version of WAM2Layers is installed
pip show wam2layers

# See what WAM2Layers can do for you
wam2layers --help
```

Now you are ready to proceed to the [user guide](./howtouse).
