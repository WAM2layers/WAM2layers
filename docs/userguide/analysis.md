
(analysis)=
# Analysis and visualization

It can be useful to quickly make some standard plots of the input and output.
To this end, WAM2layers includes a few visualization scripts.

```
# See what's available
wam2layers visualize --help

# Visualize the pre-processed data
wam2layers visualize input config-file.yaml

# Visualize the output data after tracking
wam2layers visualize output config-file.yaml
```

The configuration file is used to find the path to the data.

Of course, you are free to analyse and plot the data completely to your own
liking. These utilities are just meant to quickly inspect the data. If you find
yourself constantly making another type of plot, and you think it might be
useful to have this as a standard plot in WAM2layers, please [reach out to us
through GitHub](https://github.com/waM2layers/waM2layers/issues/new).
