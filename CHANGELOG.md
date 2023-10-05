# Changelog
All notable changes to this project will be documented in this file.

## Unreleased

### Added

- Add tests for preprocess and backtrack workflow [#265](https://github.com/WAM2layers/WAM2layers/pull/265).

### Removed

- Remove example test floodcase from CI

### Fixed

- The "analyse output" command now excludes the first timestep when it looks for output files, as it is usually not present in backtrack output ([#266](https://github.com/WAM2layers/WAM2layers/pull/266)).

### Changed

- There's now a dedicated function to compute the advection term. This makes it easier to apply different solvers. Padding is added to the edges of the domain, and boundary losses are added to e-track at the edges of the domain ([#266](https://github.com/WAM2layers/WAM2layers/pull/266)).

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
