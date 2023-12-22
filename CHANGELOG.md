# Changelog

All notable changes to this project will be documented in this file.

## Unreleased

### Added

- Add forward tracking possibility including small bugfix in backtrack ([#289](https://github.com/WAM2layers/WAM2layers/pull/289)).
- Add tests for preprocess and backtrack workflow ([#265](https://github.com/WAM2layers/WAM2layers/pull/265)).
- Add tests for visualize workflow ([#271](https://github.com/WAM2layers/WAM2layers/pull/271)).
- Add `pre-commit` to check and fix code with `black` and `isort` formatter ([#273](https://github.com/WAM2layers/WAM2layers/pull/273)).
- Add additional documentation and doctest for vertical advection terms ([#274](https://github.com/WAM2layers/WAM2layers/pull/274)).
- Copy config yaml file to output path ([#249](https://github.com/WAM2layers/WAM2layers/pull/249)).
- Add forward tracking ([#289](https://github.com/WAM2layers/WAM2layers/pull/289))
- Support for time-dependent tagging region ([#297](https://github.com/WAM2layers/WAM2layers/pull/297))
- tagging_region can be specified as a bounding box ([#297](https://github.com/WAM2layers/WAM2layers/pull/297))
- tracking_domain can now be used to subselect preprocessed data ([#297](https://github.com/WAM2layers/WAM2layers/pull/297))

### Removed

- Remove example test floodcase from CI ([#265](https://github.com/WAM2layers/WAM2layers/pull/265)).
- Chunks argument in config ([#290](https://github.com/WAM2layers/WAM2layers/pull/290))

### Fixed

- The "analyse output" command now excludes the first timestep when it looks for output files, as it is usually not present in backtrack output ([#266](https://github.com/WAM2layers/WAM2layers/pull/266)).
- A long-existing bug in the calculation of vertical fluxes, effectively changing the direction of the vertical flux ([#274](https://github.com/WAM2layers/WAM2layers/pull/274)). Also use the *new* state in the error correction term `S_1 / S_T * err_T`
- Log files are now saved to the output path, instead of the current working directory ([#291](https://github.com/WAM2layers/WAM2layers/pull/291))


### Changed

- There's now a dedicated function to compute the advection term. This makes it easier to apply different solvers. Padding is added to the edges of the domain, and boundary losses are added to e-track at the edges of the domain ([#266](https://github.com/WAM2layers/WAM2layers/pull/266)).
- Seperated the vertical advection and dispersion terms: dispersion is modelled as `kvf * abs(Fv) * dS/dz` ([#274](https://github.com/WAM2layers/WAM2layers/pull/274)).
- Splitted the backtracking module, moving functions to separate io, core, and preprocessing.shared modules ([#274](https://github.com/WAM2layers/WAM2layers/pull/274)).
- `region` renamed to `tagging_region`; `track_start_date` to `tracking_start_date`; `track_end_date` to `tracking_end_date` ([#297](https://github.com/WAM2layers/WAM2layers/pull/297))

## Release v3.0.0-beta.5 (2023-07-21)

### Added

- flowchart for users on the readme
- Improved templates of discussion forum
- source-region designer

### Removed

- Old preprocessing script from EC-Earth

### Fixed

- pydantic issue
- python versions >3.9 accepted
- last time step is outputted (but with name of 00.10)
- Output filenames now also include hours and minutes


## Release v3.0.0-beta.4 (2023-04-21)

### Added

- Documentation for developers ([#172](https://github.com/WAM2layers/WAM2layers/pull/172))
- Check that dates are in correct order ([#204](https://github.com/WAM2layers/WAM2layers/pull/204))

### Removed

- EC-Earth preprocessing script ([#195](https://github.com/WAM2layers/WAM2layers/pull/195))

### Fixed

- Datetime fields in example config file ([#194](https://github.com/WAM2layers/WAM2layers/pull/194))
- Included test of config file that start time is earlier then end time ([#204](https://github.com/WAM2layers/WAM2layers/pull/204))

## Release v3.0.0-beta.3 (2022-12-02)

### Fixed (patch version, bugs)
- n/a

### Added (minor version only)
Some highlights:
- Chunking using xarray the preprocessed data
- Visualizing output data
- Allowing global input data to read
- Timestepping of backtracking per time step (instead of double loops per day and per time within day)
- From double to single precision

### Changed and removed (major version only)
- n/a
