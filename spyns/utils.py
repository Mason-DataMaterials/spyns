# -*- coding: utf-8 -*-
"""Module with utility functions for simulations."""

__author__ = "James Glasbrenner"
__maintainer__ = "James Glasbrenner"
__email__ = "jglasbr2@gmu.edu"
__date__ = "September 8, 2017"

import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def debug_logmaker():
    """Set up logging features for debugging programs."""
    log_directory = Path.cwd()
    filename = log_directory / "debug.log"

    try:
        log_directory.mkdir()
    except FileNotFoundError:
        log_directory.mkdir(parents=True)
    except FileExistsError:
        pass

    return logging.handlers.RotatingFileHandler(
        "{0}".format(filename),
        maxBytes=31457280,
        backupCount=20,
        encoding='utf8',
    )
