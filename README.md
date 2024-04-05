[![DOI](https://zenodo.org/badge/471007521.svg)](https://zenodo.org/badge/latestdoi/471007521)
[![Documentation Status](https://readthedocs.org/projects/wam2layers/badge/?version=latest)](https://wam2layers.readthedocs.io/en/latest/?badge=latest)
[![PyPI](https://img.shields.io/pypi/v/wam2layers)](https://pypi.org/project/wam2layers/)

# Atmospheric moisture tracking

https://user-images.githubusercontent.com/17080502/233834536-a82ca96d-e303-4592-a588-472097ebe6c5.mp4

This repository contains the source code of the WAM2layers moisture tracking
code. It can be used to determine where
precipitation originally evaporated (backtracking), or where evaporated moisture
eventually ends up (forward tracking).

The animation above illustrates the concept of backtracking: you first see the water
content and fluxes move forward in time (left panel). Midway through, the animation
reverses and the moisture from the "source region" is tracked backward in time (right panel).


This code is currently developed by a core team:
Ruud van der Ent (Delft University of Technology)
Imme Benedict (Wageningen University)
Chris Weijenborg (Wageningen University)
Peter Kalverla (Netherlands eScienceCenter)
Bart Schilperoort (Netherlands eScience Center)

We are proudly part of [open-sustainable-technology](https://github.com/protontypes/open-sustainable-technology).

# How to use
See the [documentation](https://wam2layers.readthedocs.io/en/latest/) for a more detailed description.
Are you relatively new to setting up Python environments, command line interfaces etc.? Try this [video](https://youtu.be/QudlILlZnOg)

Still questions? Follow the flowchart below.
![FlowChart GitHub](https://github.com/WAM2layers/WAM2layers/assets/123247866/f5cbcf8f-a45f-4e73-b304-00956b4e2ee5)


# Other versions

This is the official codebase for the WAM-2layers moisture tracking model as of
18/03/2022, but there are still several other (older) versions around:

- [Original Python code for ERA-Interim by Ruud van der Ent](https://github.com/ruudvdent/WAM2layersPython)
- [Adapted version for EC-Earth by Imme Benedict](https://github.com/Imme1992/moisture_tracking_mississippi)
- [Adapted version for MERRA2 by Pat Keys](https://github.com/pkeys/WAM2layersPythonMerra2)
- [Adapted version for ERA5 pressure levels by Mingzhong Xiao](https://zenodo.org/record/4796962#.Y25d1-TMIVA)
- [Adapted version for ERA5 by Theo Carr](https://github.com/ktcarr/WAM2layers_ERA5)

# Reuse and acknowledgement

We are actively developing the code at the moment, so it may be subject to
change. We aim for backward compatability from v3.0.0 onward. We encourage anyone who is interested in re-using the code to get in
touch on the discussion pages. We may be able to help and you may be able to help us getting rid off bugs.

If you use the code for a publication, please cite it using the [DOI of the
appropriate release](https://doi.org/10.5281/zenodo.7010594) and the
following paper:
[Contrasting roles of interception and transpiration in the
hydrological cycle - Part 2: Moisture
recycling](https://doi.org/10.5194/esd-5-471-2014)
