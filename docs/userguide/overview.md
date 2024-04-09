(cli)=
# Using WAM2layers

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
