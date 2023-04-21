# Changelog
All notable changes to this project will be documented in this file.

## Unreleased

### Added

- Documentation for developers ([#172](https://github.com/WAM2layers/WAM2layers/pull/172))
- Check that dates are in correct order ([#204](https://github.com/WAM2layers/WAM2layers/pull/204))

### Removed

- EC-Earth preprocessing script ([#195](https://github.com/WAM2layers/WAM2layers/pull/195))

### Fixed

- Datetime fields in example config file ([#194](https://github.com/WAM2layers/WAM2layers/pull/194))


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
