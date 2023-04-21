import pytest

from wam2layers.config import Config

EXAMPLE_CONFIG = "example_config.yaml"


def test_load_example_config():
    Config.from_yaml(EXAMPLE_CONFIG)


def test_config_fails_for_wrong_time():
    config = Config.from_yaml(EXAMPLE_CONFIG)

    with pytest.raises(ValueError) as excinfo:
        config.preprocess_end_date = config.preprocess_start_date

    assert "End date should be later than start date" in str(excinfo.value)
