# Atmospheric moisture tracking

This repository contains the source code of the WAM2layers moisture tracking
code originally developed by Ruud van der Ent. It can be used to determine where
precipitation originally evaporated (backtracking), or where evaporated moisture
eventually ends up (forward tracking).

See the [documentation](http://wam2layers.rtfd.io/) for a more detailed description.

# How to use

```
# Clone the source code repository
git clone git@github.com:WAM2layers/WAM2layers.git

# Enter the source code directory
cd WAM2layers/

# Install dependencies and activate environment
conda env create -f environment.yaml
conda activate wam2layers

# Download sample data
# TODO

# Run the example backtrack script
python backtrack.py cases/FloodCase_202107.yaml
```

## Other versions

This is the official codebase for the WAM-2layers moisture tracking model as of
18/03/2022, but there are still several other versions around:

- [Original Python code by Ruud van der Ent](https://github.com/ruudvdent/WAM2layersPython)
- [Adapted version by Imme Benedict](https://github.com/Imme1992/moisture_tracking_mississippi)

## Reuse and acknowledgement
To be completed.

We are actively developing the code at the moment, so it may be subject to
change. We encourage anyone who is interested in re-using the code to get in
touch. We may be able to help.

If you use the code for a publication, please cite it using its DOI (TODO)
and/or the following paper: [Contrasting roles of interception and transpiration
in the hydrological cycle - Part 2: Moisture
recycling](https://doi.org/10.5194/esd-5-471-2014)
