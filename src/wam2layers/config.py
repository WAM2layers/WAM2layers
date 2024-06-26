import logging
from datetime import datetime
from pathlib import Path
from typing import List, Literal, NamedTuple, Optional, Union

import yaml
from pydantic import (
    AfterValidator,
    BaseModel,
    ConfigDict,
    Field,
    FilePath,
    field_validator,
    model_validator,
)
from typing_extensions import Annotated

logger = logging.getLogger(__name__)


class BoundingBox(NamedTuple):
    west: Annotated[float, Field(ge=-180, le=180)]
    south: Annotated[float, Field(ge=-80, le=80)]
    east: Annotated[float, Field(ge=-180, le=180)]
    north: Annotated[float, Field(ge=-80, le=80)]


def validate_region(region):
    """Check if region is existing path or valid bbox."""
    if isinstance(region, Path):
        return region.absolute()

    assert -180 <= region.west <= 180, "longitude should be between -180 and 180"
    assert -180 <= region.east <= 180, "longitude should be between -180 and 180"
    assert -80 <= region.north <= 80, "latitude should be between -80 and 80"
    assert -80 <= region.south <= 80, "latitude should be between -80 and 80"
    assert region.south < region.north, "south should be smaller than north"
    if region.west > region.east:
        logger.info("west > east, coordinates will be rolled around meridian")
    return region


class Config(BaseModel):
    # Pydantic configuration
    model_config = ConfigDict(validate_assignment=True)

    # Configuration settings
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

    tracking_direction: Literal["forward", "backward"]
    """The tracking direction, either forward or backward.

    You have to specify if either forward or backward tracking should be performed.

    For example:

    .. code-block:: yaml

        tracking_direction: backward

    """

    tagging_region: Annotated[
        Union[FilePath, BoundingBox], AfterValidator(validate_region)
    ]
    """Subdomain delimiting the source/sink regions for tagged moisture.

    You can either specify a path that contains a netcdf file, or a bounding box
    of the form [west, south, east, north].

    The bounding box should be inside -180, -80, 180, 80; if west > south, the
    coordinates will be rolled to retain a continous longitude.

    The file should exist. If it has a time dimension, the nearest field will be
    used as tagging region, and the time should still be between
    tagging_start_date and tagging_end_date

    For example:

    .. code-block:: yaml

        tagging_region: /data/volume_2/era5_2021/tagging_region_global.nc
        tagging_region: [0, 50, 10, 55]
    """

    tracking_domain: Optional[
        Annotated[BoundingBox, AfterValidator(validate_region)]
    ] = None
    """Subdomain delimiting the region considered during tracking.

    This is useful when you have global pre-processed data but you don't need
    global tracking.

    You can specify a bounding box of the form [west, south, east, north].

    The bounding box should be inside -180, -80, 180, 80; if west > south, the
    coordinates will be rolled to retain a continous longitude.

    If it is set to `null`, then it will use full domain of preprocessed data.

    Note that you should always set period to False if you use a subdomain.

    For example:

    .. code-block:: yaml

        tracking_domain: [0, 50, 10, 55]

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

    tracking_start_date: datetime
    """Start date for tracking.

    Should be formatted as: `"YYYY-MM-DD[T]HH:MM"`. Start date < end date, even if
    backtracking.
    When backward tracking the tracking_start_date is not given as output date.

    For example:

    .. code-block:: yaml

        tracking_start_date: "2021-07-01T00:00"

    """

    tracking_end_date: datetime
    """Start date for tracking.

    Should be formatted as: `"YYYY-MM-DD[T]HH:MM"`. Start date < end date, even if
    backtracking.

    For example:

    .. code-block:: yaml

        tracking_end_date: "2021-07-15T23:00"
    """

    tagging_start_date: datetime
    """Start date for tagging.

    For tracking individual (e.g. heavy precipitation) events, you can set the
    start and end date to something different than the total tracking start and
    end date, you can also indicate the hours that you want to track. The
    start date is included.

    Should be formatted as: `"YYYY-MM-DD[T]HH:MM"`. Start date < end date, even if
    backtracking.

    For example:

    .. code-block:: yaml

        tagging_start_date: "2021-07-13T00:00"

    """

    tagging_end_date: datetime
    """End date for tagging.

    For tracking individual (e.g. heavy precipitation) events, you can set the
    start and end date to something different than the total tracking start and
    end date, you can also indicate the hours that you want to track. The
    end date is included.

    Should be formatted as: `"YYYY-MM-DD[T]HH:MM"`. Start date < end date, even if
    backtracking.

    For example:

    .. code-block:: yaml

        tagging_end_date: "2021-07-14T23:00"

    """

    input_frequency: str
    """Frequency of the raw input data.

    Used to calculated water volumes.

    For example:

    .. code-block:: yaml

        input_frequency: '1h'

    """

    timestep: int
    """Timestep in seconds with which to perform the tracking.

    The data will be interpolated during model execution. Too large timestep will
    violate CFL criterion, too small timestep will lead to excessive numerical
    diffusion and slow progress. For best performance, the input frequency
    should be divisible by the timestep.

    For example:

    .. code-block:: yaml

        timestep: 600  # timestep in seconds

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

    levels: "Union[List[int], Literal['All']]"
    """Which levels to use from the raw input data.

    A list of integers corresponding to the levels in the input data, or a
    subset thereof. Shorthand `"all"` will attempt to use all 137 ERA5 levels.

    For example:

    .. code-block:: yaml

        levels: [20,40,60,80,90,95,100,105,110,115,120,123,125,128,130,131,132,133,134,135,136,137]

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

    kvf: float
    """net to gross vertical flux multiplication parameter

    For example:

    .. code-block:: yaml

        kvf: 3

    """

    level_layer_boundary: int = 111
    """Which level to use for the layer boundary.

    Defaults to layer number 111 (ERA5).
        TODO: Check if this is a reasonable choice for boundary

    Any layer numbers greater than the specified one will be
    included in the lower layer.

    For example:

    .. code-block:: yaml

        level_layer_boundary: 111

    """

    pressure_boundary_factor: float = 0.72878581
    """Pressure level boundary multiplication factor.

    The layer boundary is placed at:
        A * P_surf + B
    Where `P_surf` is the air pressure at the surface (at this location and time),
        `A` is the `pressure_boundary_factor`,
        and `B` is the `pressure_boundary_offset`.

    Any pressure levels above this point will end up in the upper layer.
    The others in the lower layer

    For example:

    .. code-block:: yaml

        pressure_boundary_factor: 0.7

    """

    pressure_boundary_offset: float = 7438.803223
    """Pressure level boundary offset.

    The layer boundary is placed at:
        A * P_surf + B
    Where `P_surf` is the air pressure at the surface (at this location and time),
        `A` is the `pressure_boundary_factor`,
        and `B` is the `pressure_boundary_offset`.

    Any pressure levels above this point will end up in the upper layer.
    The others in the lower layer

    For example:

    .. code-block:: yaml

        pressure_boundary_offset: 7440.0

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

    @field_validator("preprocessed_data_folder", "output_folder")
    @classmethod
    def _expanduser(cls, path):
        """Resolve ~ in input paths."""
        return path.absolute()

    @field_validator("preprocessed_data_folder", "output_folder")
    @classmethod
    def _make_dir(cls, path):
        """Create output dirs if they don't exist yet."""
        if not path.exists():
            logger.info(f"Creating output folder {path}")
            path.mkdir(parents=True)
        return path

    @model_validator(mode="after")
    def check_date_order(self):
        if self.tracking_start_date > self.tracking_end_date:
            raise ValueError(
                "tracking_end_date should be later than tracking_start_date"
            )
        if self.tagging_start_date > self.tagging_end_date:
            raise ValueError("tagging_end_date should be later than tagging_start_date")
        if self.preprocess_start_date > self.preprocess_end_date:
            raise ValueError(
                "preprocess_end_date should be later than preprocess_start_date"
            )

        if self.tracking_domain is not None and self.periodic_boundary:
            logger.warning(
                "Periodic boundary set to true while using a subdomain. Are you sure?"
            )

        return self

    def to_file(self, fname: Union[str, Path]) -> None:
        """Export the configuration to a file.

        Note that any comments and formatting from an original yaml file is lost.
        """
        fpath = Path(fname)
        with fpath.open("w") as f:
            yaml.dump(self.model_dump(mode="json"), f)
