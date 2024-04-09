# Installation

This page describes how to install WAM2Layers for various use cases, from using
the model as is, to more advanced use cases for when you want to make changes to
the code. By the end of this walkthrough you should have a working installation
of WAM2Layers and understand how you ended up there.

## Requirements

To run WAM2Layers you need a computer and some kind of command line, such as the
default terminal included in Mac or Linux. On Windows we recommend using Ubuntu
through [WSL](https://learn.microsoft.com/en-us/windows/wsl/install). If you're
using existing infrastructure such as an managed JupyterLab environment or a
supercomputer you will also typically have access to a terminal.

You also need an existing Python environment with pip (the default installer for
Python packages; often it comes included with Python), or conda/mamba to create
a new python environment for using wam2layers.

To check if Python and pip are available on your system, you can
use `which`:

```sh
# Check if python and pip are installed
which python3
which pip
```

If this returns nothing, you don't have Python/pip. In that case, we recommend
installing Python with conda. Use `which conda` to see if conda available on
your system. If it's not,
[here](https://docs.anaconda.com/free/miniconda/index.html) is a really nice
walkthrough for installing miniconda.

(pip-install)=
## Install WAM2layers in an existing Python environment

If you do have Python and pip installed, you can install wam2layers with one of the
following commands, depending on your requirements:

```
pip install wam2layers
pip install wam2layers[viz]  # includes packages required for visualization
pip install wam2layers[dev]  # include all packages needed for working on the code
pip install wam2layers=3.0.0  # install a specific version (e.g. to reproduce previous results)
```

The available versions are listed
[here](https://pypi.org/project/wam2layers/#history).

## Install WAM2layers in a fresh conda environment

If you don't have Python/pip, you can first create a new conda environment. An
additional advantage is that some packages (e.g. cartopy) might not be available
from PyPI on all systems.

To set up a dedicated environment for your moisture tracking experiments do:

```sh
# Create a fresh conda environment with cartopy, and activate it
conda create --name wamenv -c conda-forge python=3.11 cartopy
conda activate wamenv
```

Now you can install it from PyPI [as shown above](pip-install), for example:
```sh
pip install wam2layers
```

(source-install)=
## Install WAM2layers from source
In principle, the previous steps are enough to start using WAM2Layers. In our
experience, however, it is often useful and insightful to obtain a copy of the
source code and create an editable installation. This will allow you to look at
the source code to see what is going on inside the model. You can then also make
(experimental) changes to the model and see how it affects the results.

The source code for WAM2Layers is hosted on GitHub. Creating a local copy
requires git (tip: `which git`).

```sh
# Clone the source code repository
git clone https://github.com/WAM2layers/WAM2layers.git

# Enter the source code directory
cd WAM2layers/
```

Beware that the code on GitHub might include recent updates that have
not yet been included in the latest official release. This can be an advantage,
but also a downside as it may be a bit more unstable.

Local code can also be installed with pip, but instead of pointing to PyPI as
above, point to the current working directory (denoted by `.`). Additionally,
add the `--editable` flag, and `[complete]` to install all optional packages.

```
# Make an editable installation of the source code in this directory (.)
pip install --editable .[complete]
```

Now, if you edit the source code of WAM2layers, the model will run with your
updated code.

If you want to contribute your changes, raise an issue, request new features, et
cetera (highly encouraged!), check out the [contribution guidelines](./develop)
chapter of the documentation.

## Verify your installation

To confirm that WAM2layers is installed correctly, and to see which version is
installed, you can use the following commands:

```
# See if WAM2layers is installed on your system
which wam2layers

# See which version of WAM2Layers is installed
wam2layers --version
```
