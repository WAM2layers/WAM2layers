from pydantic import BaseModel, FilePath, validator, root_validator
import yaml
from datetime import datetime
from pathlib import Path
from typing import Literal, Union

import yaml
from pydantic import BaseModel, FilePath, validator


class Config(BaseModel):
    filename_template: str
    """The filename pattern of the raw input data.

    Used to find the input files during preprocessing. The pattern will be
    interpreted during execution of the model to find the input data for each
    date and variable.

    For example, the following pattern:

    .. code-block:: yaml

        filename_template: /ERA5data/{year}/{month:02}/ERA5_{year}-{month:02d}-{day:02d}{levtype}_{variable}.nc

    will be converted to

    .. code-block:: python

        /ERA5data/2021/07/ERA5_2021-07-15_ml_u.nc

    for date 2022-07-15, variable u, and levtype "_ml" (note the underscore).
    """

    preprocessed_data_folder: Path
    """Location where the pre-processed data should be stored.

    If it does not exist, it will be created during pre-processing.

    For example:

    .. code-block:: yaml

        preprocessed_data_folder: ~/floodcase_202107/preprocessed_data

    """

    region: FilePath
    """Location where the region mask is stored.

    The file should exist.

    For example:

    .. code-block:: yaml

        region: /data/volume_2/era5_2021/source_region_global.nc

    """

    output_folder: Path
    """Location where output of tracking and analysis should be written.

    For example:

    .. code-block:: yaml

        output_folder: ~/floodcase_202107/output_data

    """

    preprocess_start_date: datetime
    """Start date for preprocessing.

    Should be formatted as: `"YYYY-MM-DD[T]HH:MM"`. Start date < end date.
    The preprocess_start_date is included in the preprocessing.

    For example:

    .. code-block:: yaml

        preprocess_start_date: "2021-07-01T00:00"

    """

    preprocess_end_date: datetime
    """End date for preprocessing.

    Should be formatted as: `"YYYY-MM-DD[T]HH:MM"`. Start date < end date.
    The preprocess_end_date is included in the preprocessing.

    For example:

    .. code-block:: yaml

        preprocess_end_date: "2021-07-15T23:00"
    """

    track_start_date: datetime
    """Start date for tracking.

    Should be formatted as: `"YYYY-MM-DD[T]HH:MM"`. Start date < end date, even if
    backtracking.
    When backward tracking the track_start_date is not given as output date.

    For example:

    .. code-block:: yaml

        track_start_date: "2021-07-01T00:00"

    """

    track_end_date: datetime
    """Start date for tracking.

    Should be formatted as: `"YYYY-MM-DD[T]HH:MM"`. Start date < end date, even if
    backtracking.

    For example:

    .. code-block:: yaml

        track_end_date: "2021-07-15T23:00"
    """

    event_start_date: datetime
    """Start date for event.

    For tracking individual (e.g. heavy precipitation) events, you can set the
    start and end date of the event to something different than the total
    tracking start and end date, you can also indicate the hours that you want to track.
    The event_start_date is included in the event tracking.

    Should be formatted as: `"YYYY-MM-DD[T]HH:MM"`. Start date < end date, even if
    backtracking.

    For example:

    .. code-block:: yaml

        event_start_date: "2021-07-13T00:00"

    """

    event_end_date: datetime
    """Start date for event.

    For tracking individual (e.g. heavy precipitation) events, you can set the
    start and end date of the event to something different than the total
    tracking start and end date, you can also indicate the hours that you want to track.
    The event_end_date is included in the event tracking.

    Should be formatted as: `"YYYY-MM-DD[T]HH:MM"`. Start date < end date, even if
    backtracking.

    For example:

    .. code-block:: yaml

        event_end_date: "2021-07-14T23:00"

    """

    input_frequency: str
    """Frequency of the raw input data.

    Used to calculated water volumes.

    For example:

    .. code-block:: yaml

        input_frequency: '1h'

    """

    target_frequency: str
    """Frequence at which to perform the tracking.

    The data will be interpolated during model execution. Too low frequency will
    violate CFL criterion, too high frequency will lead to excessive numerical
    diffusion and slow progress. For best performance, the input frequency
    should be divisible by the target frequency.

    For example:

    .. code-block:: yaml

        target_frequency: '10min'

    """

    output_frequency: str
    """Frequency at which to write output to file.

    For example, for daily output files:

    .. code-block:: yaml

        output_frequency: '1d'

    """

    level_type: Literal["model_levels", "pressure_levels"]
    """Type of vertical levels in the raw input data.

    Can be either `model_levels` or `pressure_levels`.

    For example:

    .. code-block:: yaml

        level_type: model_levels

    """

    levels: "Union[list[int], Literal['All']]"
    """Which levels to use from the raw input data.

    A list of integers corresponding to the levels in the input data, or a
    subset thereof. Shorthand `"all"` will attempt to use all 137 ERA5 levels.

    For example:

    .. code-block:: yaml

        levels: [20,40,60,80,90,95,100,105,110,115,120,123,125,128,130,131,132,133,134,135,136,137]

    """

    log_level: Literal["debug", "info", "warning", "error", "critical"]
    """Verbosity of the output messages.

    For example:

    .. code-block:: yaml

        log_level: info
    """

    restart: bool
    """Whether to restart from previous run.

    If set to `true`, this will attempt to read the output from a previous model
    run and continue from there. The output from the previous timestep must be
    available for this to work.

    For example:

    .. code-block:: yaml

        restart: False

    """

    periodic_boundary: bool
    """Whether to use period boundaries in the zonal direction.

    This should be used when working with global datasets.

    For example:

    .. code-block:: yaml

        periodic_boundary: true
    """

    kvf: int
    """net to gross vertical flux multiplication parameter

    For example:

    .. code-block:: yaml

        kvf: 3

    """

    chunks: "Union[None, dict[str, int]]"
    """Whether to use dask.

    Using dask can help process large datasets without memory issues, but its
    performance is quite sensitive to the chunk configuration.

    Some examples:

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

    @validator("preprocessed_data_folder", "region", "output_folder")
    def _expanduser(cls, path):
        """Resolve ~ in input paths."""
        return path.expanduser()

    @validator("preprocessed_data_folder", "output_folder")
    def _make_dir(cls, path):
        """Create output dirs if they don't exist yet."""
        if not path.exists():
            print(f"Creating output folder {path}")
            path.mkdir(parents=True)
        return path

    @root_validator
    def check_date_order(cls, values):
        for date_field in ['track', 'event', 'preprocess']:
            start = values.get(f'{date_field}_start_date')
            end = values.get(f'{date_field}_end_date')

            if not start < end:
                raise ValueError("End date should be later than start date.")

        return values

    class Config:
        """This is the config of the **Pydantic** class, not the wam2layers config.

        Note: this will change in v2
        https://docs.pydantic.dev/blog/pydantic-v2-alpha/#changes-to-config
        """
        # Also force validation when new value is set on existing config
        validate_assignment = True

