(interactive)=
## Using WAM2layers in an interactive session

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
