from pydantic import BaseModel, FilePath
import yaml
from datetime import date
from typing import Literal, Union
from pathlib import Path


class Config(BaseModel):
    # Settings needed for the data preprocessing
    filename_template: str
    """The filename pattern of the raw input data.

    Used to find the input files during preprocessing. The pattern will be
    interpreted during execution of the model to find the input data for each
    date and variable. For example, the following pattern:

    `/ERA5data/{year}/{month:02}/ERA5_{year}-{month:02d}-{day:02d}{levtype}_{variable}.nc`

    will be converted to

    `/ERA5data/2021/07/ERA5_2021-07-15_ml_u.nc`

    for date `2022-07-15`, variable `u` and levtype `_ml` (note the underscore).
    """

    preprocess_start_date: date
    """Start date for preprocessing.

    Should be formatted as: `"YYYY-MM-DD"`. Start date < end date.
    """

    preprocess_end_date: date
    """End date for preprocessing.

    Should be formatted as: `"YYYY-MM-DD"`. Start date < end date.
    """

    level_type: Literal["model_levels", "pressure_levels"]
    """Type of vertical levels in the raw input data.

    Can be either `model_levels` or `pressure_levels`.
    """

    levels: Union[list[int], Literal["All"]]
    """Which levels to use from the raw input data.

    A list of integers corresponding to the levels in the input data, or a
    subset thereof. Shorthand `"all"` will attempt to use all 137 ERA5 levels.
    """

    # Settings shared by preprocessing and backtracking
    preprocessed_data_folder: Path
    """Location where the pre-processed data should be stored.

    For example: `~/floodcase_202107/preprocessed_data`

    If it does not exist, it will be created during pre-processing.
    """

    # Settings needed to define the tracking region
    region: FilePath
    """Location where the region mask is stored.

    For example: `/data/volume_2/era5_2021/source_region_global.nc`
    The file should exist.
    """

    track_start_date: date
    """Start date for tracking.

    Should be formatted as: `"YYYY-MM-DD"`. Start date < end date, even if
    backtracking.
    """

    track_end_date: date
    """Start date for tracking.

    Should be formatted as: `"YYYY-MM-DD"`. Start date < end date, even if
    backtracking.
    """

    # Settings needed for the tracking run
    input_frequency: str
    """Frequency of the raw input data.

    Used to calculated water volumes. For example: `'1h'`
    """

    target_frequency: str
    """Frequence at which to perform the tracking.

    For example: `'15m'`.

    The data will be interpolated during model execution. Too low frequency will
    violate CFL criterion, too high frequency will lead to excessive numerical
    diffusion and slow progress. For best performance, the input frequency
    should be divisible by the target frequency.
    """

    output_frequency: str
    """Frequency at which to write output to file.

    For example: `'1d'` for daily output files.
    """

    periodic_boundary: bool
    """Whether to use period boundaries in the zonal direction.

    For example: `true`
    """

    output_folder: Path
    """Location where output of tracking and analysis should be written.

    For example:  `~/floodcase_202107/output_data`
    """

    restart: bool
    """Whether to restart from previous run.

    If set to `true`, this will attempt to read the output from a previous model
    run and continue from there. The output from the previous timestep must be
    available for this to work.
    """

    kvf: int
    """Stability correction parameter.

    For example: `3`
    """

    timetracking: bool
    """Whether to also track residence time of parcels.

    Currently not implemented.
    """

    distancetracking: bool
    """Whether to also track distance traveled by parcels.

    Currently not implemented.
    """

    log_level: Literal["debug", "info", "warning", "error", "critical"]
    """Verbosity of the output messages.

    Set to `debug` to get the most information printed during model executing.
    """

    chunks: Union[None, dict[str, int]]
    """Whether to use dask.

    Using dask can help process large datasets without memory issues, but its
    performance is quite sensitive to the chunk configuration.

    .. highlight:: yaml
    .. code-block:: yaml

        # don't use dask:
        chunks: null

        # process one time step at a time
        chunks: {time: 1}

        # set chunking for all individual dims
        chunks:
            level: -1
            time: 1  # don't set time to auto, as it will be different for surface and 3d vars
            latitude: auto
            longitude: auto

    """

    event_start_date: date
    """Start date for event.

    For tracking individual (e.g. heavy precipitation) events, you can set the
    start and end date of the event to something different than the total
    tracking start and end date.

    Should be formatted as: `"YYYY-MM-DD"`. Start date < end date, even if
    backtracking.
    """

    event_end_date: date
    """Start date for event.

    For tracking individual (e.g. heavy precipitation) events, you can set the
    start and end date of the event to something different than the total
    tracking start and end date.

    Should be formatted as: `"YYYY-MM-DD"`. Start date < end date, even if
    backtracking.
    """

    @classmethod
    def from_yaml(cls, config_file):
        """Read settings from a configuration.yaml file.

        For example:

        .. highlight:: python
        .. code-block:: python

            from wam2layers.config import Config
            config = Config.from_yaml('../../cases/floodcase_2021.yaml')

        """
        with open(config_file) as f:
            settings = yaml.safe_load(f)

        return cls(**settings)
