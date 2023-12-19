"""Utils for logging."""
import logging
from datetime import datetime
from pathlib import Path


def setup_logging(log_path):
    """Configure logging behaviour.

    Messages with logger.INFO (and higher) are written to stdout
    Messages with logger.DEBUG (and higher) are written to wam2layers.log
    """
    # https://stackoverflow.com/a/58828499
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")

    file_handler = logging.FileHandler(
        Path(log_path, f"wam2layers_{timestamp}.log"), mode="w"
    )
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    logging.basicConfig(
        level=logging.INFO,
        format="%(message)s",
        handlers=[file_handler, stream_handler],
    )
