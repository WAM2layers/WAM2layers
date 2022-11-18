# How to use

## Quickstart
We recommend to create an editable installation in a dedicated conda environment (assuming you have an installation of Anaconda/Miniconda already):

```
# Create a dedicated conda environment
conda create --name wam2env -c conda-forge python=3.9 jupyterlab cartopy matplotlib cmocean
conda activate wam2env

# Clone the source code repository
git clone git@github.com:WAM2layers/WAM2layers.git
cd WAM2layers

# Install the package from source
pip install --editable .

# Download example data and update the paths in the configuration file
https://rebrand.ly/wam2layers-example-data

# Run example case
wam2layers backtrack floodcase_202107.yaml
```





