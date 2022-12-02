from pydantic import BaseModel, FilePath
import yaml
from datetime import date
from typing import Literal
from pathlib import Path


class Config(BaseModel):
    # Settings needed for the data preprocessing
    filename_template: str
    preprocess_start_date: date
    preprocess_end_date: date
    level_type: Literal["model_levels", "pressure_levels"]
    levels: list[int]

    # Settings shared by preprocessing and backtracking
    preprocessed_data_folder: Path

    # Settings needed to define the tracking region
    region: FilePath
    track_start_date: date
    track_end_date: date

    # Settings needed for the tracking run
    input_frequency: str
    target_frequency: str
    output_frequency: str
    periodic_boundary: bool
    output_folder: Path
    restart: bool
    kvf: int
    timetracking: bool
    distancetracking: bool
    log_level: Literal["debug", "info", "warning", "error", "critical"]
    chunks: dict[str, int]

    event_start_date: date
    event_end_date: date

    @classmethod
    def from_yaml(cls, config_file):
        with open(config_file) as f:
            settings = yaml.safe_load(f)

        return cls(**settings)
