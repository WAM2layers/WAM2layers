(cli)=
# Using WAM2layers

## Through the command line

The primary way of running WAM2layers is through the command line, much
like pip, conda, git, cdo, et cetera. It looks like this:

```sh
wam2layers track experiment-config-file.yaml
```

To see all commands that are available from wam2layers, you can use

```sh
wam2layers --help
```

Most commands require a config file. This file contains the settings for your
experiment. More information is available under [](config).

The command line interface is very useful when you want to do multiple
experiments with different settings. You only ever edit the config file, not the
source code. It is also useful for running the model on an HPC system. Moreover,
the command `wam2layers` can be used from anywhere on your system.

```{Important}
A possible disadvantage of the command line interface is that it obscures what
is happening under the hood. Note that it is also possible to download the source
code, make modifications, and run your custom version of WAM2layers instead.

Once you get used to working with WAM2layers, we strongly encourage you to
[create such a source installation](source-install) as well.
```


(interactive)=
## In an interactive session

While WAM2layers was designed to be run from the command line, it is possible to
run it in an interactive session such as IPython or Jupyter Lab. You can import
WAM2layers functionality and build your own scripts. For example:

```python
import wam2layers
config_file = "path_to/example_data/floodcase_2021.yaml"
wam2layers.tracking.backtrack.run_experiment(config_file)
```

You can also import specific functions. However, you should be conscious that
the code was not designed to be used in this way, so it might not be the most
intuitive, and you might run into unexpected behaviour. Moreover, in maintaining
the model, we cannot guarantee backward compatibility for this type of use. So
use it at your own risk.

```{admonition} Idea
Perhaps, at some point, it would be nice to build a more user-friendly Python
API, for example following the ["basic model
interface"](https://bmi.readthedocs.io/en/stable/):

```python
from wam2layers import BacktrackModel

config_file = "path_to/example_data/floodcase_2021.yaml"

model = BacktrackModel()
model.initialize(config_file)

while model.time < model.end_time:
    model.update()

model.get_value('s_upper')
```
