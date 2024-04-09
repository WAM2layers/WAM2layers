# Case configuration

WAM2layers uses case configuration files to store the settings for an
experiment. That makes it possible to run various experiments without changing
the model code. The configuration files are written in
[yaml](https://yaml.org/spec/1.2.2/#chapter-1-introduction-to-yaml) format. When
you run WAM2layers, these settings are loaded into a `Config` object. The
options in your yaml file should correspond to the attributes of the `Config`
class listed below.

<!-- TODO: update links -->
Configuration files for the two example cases are available under []().


<!-- This generates automatic documentation based on docstrings in src/wam2layers/config.py -->
```{eval-rst}
.. autoclass:: wam2layers.config.Config
    :members:
    :undoc-members:
```
