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
wam2layers preprocess floodcase_202107.yaml # will not work unless you also have data for the pre-processing 
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