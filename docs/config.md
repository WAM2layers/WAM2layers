## Case configuration

WAM2layers uses case configuration files to store the settings for an
experiment. That makes it possible to run various experiments without changing
the model code.

The configuration files are written in [yaml]() format and can/should include the following settings:

**Settings needed for the data preprocessing (as currently for era5)**
- `preprocess_start_date`: formatted as 'YYYYMMDD', e.g. '20210701'
- `preprocess_end_date`: formatted as 'YYYYMMDD', e.g. '20210716'
- `level_type`: either "pressure_levels" or "model_levels"
- `levels`: "all" or a list of integers with the desired (model or pressure)
  levels.
- `input_folder`: path where raw era5 input data can be found, e.g.
  /home/peter/WAM2layers/era5_2021
- `filename_prefix`: Fixed part of filename. This function will infer the variable
  name and add _ml for model level data. E.g. with prefix = "FloodCase_202107"
  this function will be able to find FloodCase_202107_ml_u.nc or
  FloodCase_202107_u.nc and FloodCase_202107_sp.nc

**Settings shared by preprocessing and backtracking**
- `preprocessed_data_folder`: path where preprocessed data should be stored. This
  directory will be created if it does not exist. E.g.
  /home/peter/WAM2layers/preprocessed_data_2021

**Settings needed to define the tracking region (in space and time)**
- `region`: path to source region. Must have a variable named "source_region",
  which must have the same extent in space as preprocessed data and values
  between 0 and 1. E.g. /home/peter/WAM2layers/era5_2021/source_region.nc
- `track_start_date`: formatted as 'YYYYMMDD', e.g. '20210701'
- `track_end_date`: formatted as 'YYYYMMDD', e.g. '20210716'

**Settings needed for the tracking run**
- `target_frequency`: WAM2layers interpolates the data internally to avoid
  numerical instabilities. Use a frequency string, e.g. '15min'. See
  https://stackoverflow.com/a/35339226 for options
- `periodic_boundary`: true or false. Set to true if input data goes from 180W to
  180E
- `output_folder`: path to output folder, e.g.
  /home/peter/WAM2layers/output_data_2021
- `restart`: true or false. Whether to load tracked water from previous run. If
  false, start from zero tracked water
- `kvf`: Vertical transport parameter for gross vertical transport between the
  layers during the tracking: "actual exchange = Kvf * F_vertical + F_vertical"
  in one direction and "-1 * (Kvf * F_vertical)" in opposite direction. A good
  default choice is kvf=3.
- `timetracking`: true or false. Currently nothing is done with this setting
- `distancetracking`: true or false. Currently nothing is done with this setting
- `log_level`: "high" or "low". Setting it too high will print more info/warning messages.

**Settings that control the period in which water is tracked**
- `event_start_date`: formatted as 'YYYYMMDD', e.g. '20210713'
- `event_end_date`: formatted as 'YYYYMMDD', e.g. '20210715'
