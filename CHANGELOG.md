# Changelog

All notable changes to this project will be documented in this file.

## Unreleased

### Added

- Windows instruction to the installation guide ([#392](https://github.com/WAM2layers/WAM2layers/pull/392)).
- support for preprocessing CMIP data ([#401](https://github.com/WAM2layers/WAM2layers/pull/401)).
- experimental support for parallel preprocessing ([#401](https://github.com/WAM2layers/WAM2layers/pull/401)).
- config file documentation expanded ([#398](https://github.com/WAM2layers/WAM2layers/pull/398)).
- support for preprocessing ARCO-ERA5 data ([#401](https://github.com/WAM2layers/WAM2layers/pull/401)).
- longitude shifter to preprocessing that ensures all longitude values are in the range (-180, 180) ([#401](https://github.com/WAM2layers/WAM2layers/pull/401)).

### Changed
- input data to the tracking code is now assumed to have longitude values between -180 and 180 degrees ([#401](https://github.com/WAM2layers/WAM2layers/pull/401)).

### Fixed

- Patch for bug in meridional advection not accounting for decreasing grid cell size towards poles ([#320](https://github.com/WAM2layers/WAM2layers/pull/320)).
- In simulations with periodic boundaries, stagger the flux to the exterior edges as well ([#426](https://github.com/WAM2layers/WAM2layers/pull/426)).
- support for ERA5 data from new CDS ([#429](https://github.com/WAM2layers/WAM2layers/pull/429)).

## Release v3.1.0 (2024-06-21)

### Changed

- The user guide has been restructured [#367](https://github.com/WAM2layers/WAM2layers/pull/367)
- The example config has been updated [#371](https://github.com/WAM2layers/WAM2layers/pull/371)
- Added config options to allow for specifying the layer boundary during preprocessing [#391](https://github.com/WAM2layers/WAM2layers/pull/391).
- The preprocessing module has been restructured to allow for easier implementation of new datasets [#391](https://github.com/WAM2layers/WAM2layers/pull/391)

### Fixed

- Fixed bug in treatment of periodic boundaries ([#379](https://github.com/WAM2layers/WAM2layers/pull/379)).
- Pressure level (ERA5) data can be preprocessed again [#387](https://github.com/WAM2layers/WAM2layers/pull/387)

## Release v3.0.0 (2024-04-05)

### Added

- Updated the visualization: more generic for forward and backward tracking ([#318](https://github.com/WAM2layers/WAM2layers/pull/318)).
- Add difference plots to the regression tests (`test_workflow.py`). The plots are uploaded as "artifacts" on Github Actions for easy inspection ([#319](https://github.com/WAM2layers/WAM2layers/pull/319)).
- Added an export to file method for the configuration object (`wam2layers.config.Config`) to allow for easier testing of the command line interface ([#320](https://github.com/WAM2layers/WAM2layers/pull/320)).
- Added a regression test for the forward tracking workflow ([#320](https://github.com/WAM2layers/WAM2layers/pull/320)).
- Added support for Python 3.12 ([#333](https://github.com/WAM2layers/WAM2layers/pull/333)).
- Adopted [NEP-29](https://numpy.org/neps/nep-0029-deprecation_policy.html) as version support policy ([#333](https://github.com/WAM2layers/WAM2layers/pull/333)).
- The output files now contain many attributes for easier interpretation ([#334](https://github.com/WAM2layers/WAM2layers/pull/334)).
- Publishing of the package to the Python Package Index ([PyPI](https://pypi.org/)) is now automated with a Github Actions workflow ([#342](https://github.com/WAM2layers/WAM2layers/pull/342)).
- You can now view the version of wam2layers by running `wam2layers --version` ([#352](https://github.com/WAM2layers/WAM2layers/pull/352)).
- Added a debug flag to the CLI: `wam2layers --debug`. This will make debug level statements be printed to the terminal as well ([#354](https://github.com/WAM2layers/WAM2layers/pull/354)).

### Removed

- Support for Python 3.8 has been removed. This is in line with the new policy ([following NEP-29](https://numpy.org/neps/nep-0029-deprecation_policy.html)) ([#333](https://github.com/WAM2layers/WAM2layers/pull/333)).
- The CLI command to generate snapshots was removed. This was deemed to unstable for a V3 release ([#362](https://github.com/WAM2layers/WAM2layers/pull/362)).

### Fixed

- Fixed a bug in the profiler's `track_stability_correction` method ([#325](https://github.com/WAM2layers/WAM2layers/pull/325)).
- Log message about "lost moisture" no longer includes boundary losses. Instead,
  there is a separate log message for "boundary transport"
  ([#343](https://github.com/WAM2layers/WAM2layers/pull/343))
- Fixed a bug in the calculation of losses and gains for backtrack, where the moisture was compared to the wrong reference state ([#355](https://github.com/WAM2layers/WAM2layers/pull/355)).
- The calculated gains are now absolute fields instead of relative to the reference field ([#355](https://github.com/WAM2layers/WAM2layers/pull/355)).

### Changed

- The workflow tests use a temporary directory managed by pytest and the user's operating system, to avoid the `/tmp` folder polluting the workspace ([#320](https://github.com/WAM2layers/WAM2layers/pull/320)).
- The workflow tests now make extensive use of pytest's [fixtures](https://docs.pytest.org/en/8.0.x/explanation/fixtures.html), which will make future test writing easier ([#320](https://github.com/WAM2layers/WAM2layers/pull/320)).
- The `kvf` parameter can now be a floating point number instead of an integer ([#320](https://github.com/WAM2layers/WAM2layers/pull/320)).
- The stability correction warning will now occur only when 10% more of the grid is corrected compared to last warning, or the correction factor is 10% stronger ([#332](https://github.com/WAM2layers/WAM2layers/pull/332)).
- The output files now contain a time dimension (incl. the coordinate) ([#334](https://github.com/WAM2layers/WAM2layers/pull/334)).
- The package version is now defined in `src/wam2layers/__init__.py` ([#334](https://github.com/WAM2layers/WAM2layers/pull/334)).
- The command line interface for tracking has changed. Tracking experiments are now started by doing `wam2layers track config.yml`. The tracking direction is retrieved from the new configuration entry "tracking_direction" ([#338](https://github.com/WAM2layers/WAM2layers/pull/338)).
- The configuration parameter "target_frequency" has been renamed to "timestep". Also, under the hood, it is applied in a different location, closer to the core computations. This makes it easier to see the relation between the formula "on paper" and the implementation in the code ([#346](https://github.com/WAM2layers/WAM2layers/pull/346)).
- The tests on Github Actions now use base python, not micromamba ([#351](https://github.com/WAM2layers/WAM2layers/pull/351)).
- The units with which WAM2layers calculates internally are changed from m3 to kg/m2. This saves some unit conversions, makes it easier to see the relation between the equations and the code, and makes it easier to interpret the output and intermediate states. Visualization functions no longer need to convert back to kg/m2 ([#356](https://github.com/WAM2layers/WAM2layers/pull/356)).
- The preprocessed data now mostly follows the CF-1.6 convention for netCDF files ([#363](https://github.com/WAM2layers/WAM2layers/pull/363)).
- The output data now follows the CF-1.6 convention for netCDF files ([#363](https://github.com/WAM2layers/WAM2layers/pull/363)).

## Release v3.0.0-beta.6 (2023-12-22)

### Added

- Add information about losses to profiling and screen messages | losses and gains as output data ([#305](https://github.com/WAM2layers/WAM2layers/pull/305)).
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
- A long-existing bug in the calculation of vertical fluxes, effectively changing the direction of the vertical flux ([#274](https://github.com/WAM2layers/WAM2layers/pull/274)). Also use the _new_ state in the error correction term `S_1 / S_T * err_T`
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
