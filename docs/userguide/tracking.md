
(tracking)=
# Tracking

The core of WAM2layers are the tracking routines. The model includes forward and
backward tracking. Forward tracking takes evaporation over the tagging region as
input and generates tracked precipitation as output. Backward tracking takes
precipitation over the tagging region as input and generates tracked evaporation
as output. Time tracking, distance tracking, and moisture recycling can
be added in future updates.

Assuming you have a preprocessed dataset and prepared a configuration file for
your experiment, running WAM2layers is as simple as:

```
wam2layers track config-file.yaml
```

where `config-file.yaml` is the path to your configuration file. Among others,
this file should have settings on the tracking direction (forward/backward),
the date range for which you want to run track, and also about the location where the
preprocessed data are stored and where the output will be stored.
For more information, see [](./config) or have a look at the example config file
[here](https://github.com/WAM2layers/WAM2layers/blob/main/example_config.yaml).

```{tip}
The actual code that is executed can be found in
`src/wam2layers/tracking/backtrack.py` or `src/wam2layers/tracking/forwardtrack.py`. This script reads the configuration,
loads the preprocessed data step by step, calculates the proportion of tracked
moisture in each grid cell, and writes the output to a file. You can configure
how often output is written.

For more information on the equations being solved, see the
[theory](../theory.md) chapter. If you want to know how the command line call
ends up in this python script, see the section on [WAM2layers' command line
interface](cli)
```
