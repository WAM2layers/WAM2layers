# Case configuration

WAM2layers uses case configuration files to store the settings for an
experiment. That makes it possible to run various experiments without changing
the model code. The configuration files are written in
[yaml](https://yaml.org/spec/1.2.2/#chapter-1-introduction-to-yaml) format. When
you run WAM2layers, these settings are loaded into a `Config` object. The
options in your yaml file should correspond to the attributes of the `Config`
class listed below.

An example configuration file is available
[here](https://github.com/WAM2layers/WAM2layers/blob/main/example_config.yaml).
Alternatively, checkout the configuration files for the example cases for
[Volta](https://data.4tu.nl/datasets/bbe10a2a-39dc-4098-a69f-0f677d06ecdd) and
[Eiffel](https://data.4tu.nl/datasets/f9572240-f179-4338-9e1b-82c5598529e2).


<!-- This generates automatic documentation based on docstrings in src/wam2layers/config.py -->
```{eval-rst}
.. autoclass:: wam2layers.config.Config
    :members:
    :undoc-members:
```
