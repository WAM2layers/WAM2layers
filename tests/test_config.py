import logging
from datetime import datetime

import pytest

from wam2layers.config import BoundingBox, Config

EXAMPLE_CONFIG = "example_config.yaml"
TEST_CONFIG = "tests/test_data/config_rhine.yaml"


def test_load_example_config():
    """Verify that the example config file can be loaded"""
    Config.from_yaml(EXAMPLE_CONFIG)


def test_load_test_config():
    """Verify that the test config file can be loaded"""
    Config.from_yaml(TEST_CONFIG)


@pytest.mark.parametrize("config_file", [EXAMPLE_CONFIG, TEST_CONFIG])
def test_config_export(tmp_path, config_file):
    """Verify that a config written to a file is identical to the original one."""
    cfg = Config.from_yaml(EXAMPLE_CONFIG)
    export_path = tmp_path / "config.yaml"
    cfg.to_file(export_path)
    cfg_exported = Config.from_yaml(export_path)
    assert cfg_exported == cfg


def test_check_date_order_invalid_tracking():
    """Raise error when tracking_end_date is earlier than tracking_start_date."""
    config = Config.from_yaml(TEST_CONFIG)
    with pytest.raises(
        ValueError, match="tracking_end_date should be later than tracking_start_date"
    ):
        config.tracking_end_date = datetime(2022, 7, 31)


def test_check_date_order_invalid_tagging():
    """Raise error when tagging_end_date is earlier than tagging_start_date."""
    config = Config.from_yaml(TEST_CONFIG)
    with pytest.raises(
        ValueError, match="tagging_end_date should be later than tagging_start_date"
    ):
        config.tagging_end_date = datetime(2022, 7, 31)


def test_check_date_order_invalid_preprocess():
    """Raise error when preprocess_end_date is earlier than preprocess_start_date."""
    config = Config.from_yaml(TEST_CONFIG)
    with pytest.raises(
        ValueError,
        match="preprocess_end_date should be later than preprocess_start_date",
    ):
        config.preprocess_end_date = datetime(2022, 7, 31)


class TestValidateRegion:
    def test_path_expanded(self):
        """Check path is expanded."""
        config = Config.from_yaml(TEST_CONFIG)
        assert config.tagging_region.is_absolute()
        config.tagging_region = "tests/test_data/region_rhine.nc"
        assert config.tagging_region.is_absolute()

    def test_invalid_path(self):
        """Raise error for invalid path."""
        config = Config.from_yaml(TEST_CONFIG)
        with pytest.raises(Exception, match="does not point to a file"):
            config.tagging_region = "weird/path.exe"

    def test_valid_bbox_ok(self):
        """Check that we can set region to a bbox"""
        config = Config.from_yaml(TEST_CONFIG)
        config.tagging_region = BoundingBox(0, 50, 10, 55)
        assert config.tagging_region == BoundingBox(0, 50, 10, 55)

    def test_invalid_south(self):
        """Raise error for invalid bbox."""
        config = Config.from_yaml(TEST_CONFIG)
        with pytest.raises(Exception, match="south should be smaller than north"):
            config.tagging_region = BoundingBox(0, 60, -10, 55)

    def test_wrap_meridian(self, caplog):
        """Warn for bbox crossing the antimeridian"""
        config = Config.from_yaml(TEST_CONFIG)
        with caplog.at_level(logging.INFO):
            config.tagging_region = BoundingBox(10, 50, -10, 55)
        assert "coordinates will be rolled" in caplog.text
        assert config.tagging_region.west == 10
