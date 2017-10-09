# -*- coding: utf-8 -*-
"""Module with utility functions for simulations."""

__author__ = "James Glasbrenner"
__maintainer__ = "James Glasbrenner"
__email__ = "jglasbr2@gmu.edu"
__date__ = "September 8, 2017"

import logging
import math
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
        "{0}".format(filename), maxBytes=31457280, backupCount=20,
        encoding='utf8')


def convert_spherical_to_cartesian(rho, phi, theta):
    """Convert from spherical coordinates into cartesian coordinates.

    Parameters
    ----------
    rho : float
        3D vector length

    phi : float
        Counterclockwise angle in the xy plane from the x axis (phi = 0) to the
        xy projection of the 3D vector.

    theta : float
        The angle between the z axis and the 3D vector.

    Returns
    -------
    list
        The [x, y, z] components corresponding to the vector defined by `rho`,
        `phi`, and `theta`.

    """
    # convert phi and theta into radians
    phi = phi * math.pi / 180
    theta = theta * math.pi / 180

    # Cartesian coordinates
    x = rho * math.sin(theta) * math.cos(phi) if abs(rho) >= 1E-9 else 0
    y = rho * math.sin(theta) * math.sin(phi) if abs(rho) >= 1E-9 else 0
    z = rho * math.cos(theta) if abs(rho) >= 1E-9 else 0

    return [
        0 if abs(x) < 1E-9 else x,
        0 if abs(y) < 1E-9 else y,
        0 if abs(z) < 1E-9 else z,
    ]
