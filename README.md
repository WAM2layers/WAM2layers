[![DOI](https://zenodo.org/badge/471007521.svg)](https://zenodo.org/badge/latestdoi/471007521)
[![Documentation Status](https://readthedocs.org/projects/wam2layers/badge/?version=latest)](https://wam2layers.readthedocs.io/en/latest/?badge=latest)
[![PyPI](https://img.shields.io/pypi/v/wam2layers)](https://pypi.org/project/wam2layers/)

# Atmospheric moisture tracking

This repository contains the source code of the WAM2layers moisture tracking
code. This code is currently developed by a core team:
Ruud van der Ent (Delft University of Technology)
Imme Benedict (Wageningen University)
Chris Weijenborg (Wageningen University)
Peter Kalverla (Netherlands eScienceCenter)

It can be used to determine where
precipitation originally evaporated (backtracking), or where evaporated moisture
eventually ends up (forward tracking).

# How to use
See the [documentation](http://wam2layers.rtfd.io/) for a more detailed description.

## Other versions

This is the official codebase for the WAM-2layers moisture tracking model as of
18/03/2022, but there are still several other older versions around:

- [Original Python code for ERA-Interim by Ruud van der Ent](https://github.com/ruudvdent/WAM2layersPython)
- [Adapted version for EC-Earth by Imme Benedict](https://github.com/Imme1992/moisture_tracking_mississippi)
- [Adapted version for MERRA2 by Pat Keys](https://github.com/pkeys/WAM2layersPythonMerra2)
- [Adapted version for ERA5 pressure levels by Mingzhong Xiao](https://zenodo.org/record/4796962#.Y25d1-TMIVA)

## Reuse and acknowledgement

We are actively developing the code at the moment, so it may be subject to
change. We encourage anyone who is interested in re-using the code to get in
touch. We may be able to help.

If you use the code for a publication, please cite it using the [DOI of the
appropriate release](https://doi.org/10.5281/zenodo.7010594) and the
following paper: 
[Contrasting roles of interception and transpiration in the
hydrological cycle - Part 2: Moisture
recycling](https://doi.org/10.5194/esd-5-471-2014)

