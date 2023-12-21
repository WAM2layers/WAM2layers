from datetime import datetime

import pytest

from wam2layers.config import Config

EXAMPLE_CONFIG = "example_config.yaml"
TEST_CONFIG = "tests/test_data/test_config.yaml"


def test_load_example_config():
    """Verify that the example config file can be loaded"""
    Config.from_yaml(EXAMPLE_CONFIG)


def test_load_test_config():
    """Verify that the test config file can be loaded"""
    Config.from_yaml(TEST_CONFIG)


def test_check_date_order_invalid_tracking():
    valid_config = Config.from_yaml(TEST_CONFIG)
    # Ensure ValueError is raised when tracking_end_date is earlier than tracking_start_date
    invalid_config = valid_config.model_copy(
        update={"tracking_end_date": datetime(2022, 7, 31)}
    )
    with pytest.raises(
        ValueError, match="tracking_end_date should be later than tracking_start_date"
    ):
        invalid_config.check_date_order()


def test_check_date_order_invalid_tagging():
    valid_config = Config.from_yaml(TEST_CONFIG)
    # Ensure ValueError is raised when tagging_end_date is earlier than tagging_start_date
    invalid_config = valid_config.model_copy(
        update={"tagging_end_date": datetime(2022, 7, 31)}
    )
    with pytest.raises(
        ValueError, match="tagging_end_date should be later than tagging_start_date"
    ):
        invalid_config.check_date_order()


def test_check_date_order_invalid_preprocess():
    valid_config = Config.from_yaml(TEST_CONFIG)
    # Ensure ValueError is raised when preprocess_end_date is earlier than preprocess_start_date
    invalid_config = valid_config.model_copy(
        update={"preprocess_end_date": datetime(2022, 7, 31)}
    )
    with pytest.raises(
        ValueError,
        match="preprocess_end_date should be later than preprocess_start_date",
    ):
        invalid_config.check_date_order()
