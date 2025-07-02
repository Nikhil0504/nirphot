import logging
import os
import sys


def setup_logging(level=None):
    """Set up basic logging for the package or CLI."""
    level = level or os.getenv("LOGLEVEL", "INFO").upper()
    level = getattr(logging, level, logging.INFO)

    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter(
        "[%(asctime)s] %(levelname)s - %(name)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    handler.setFormatter(formatter)

    root_logger = logging.getLogger()
    root_logger.setLevel(level)

    if not root_logger.handlers:
        root_logger.addHandler(handler)
