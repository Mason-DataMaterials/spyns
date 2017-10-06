# -*- coding: utf-8 -*-
"""Wrapper for pymatgen Structure objects.

Copyright (c) Mason DataMaterials Group.
Distributed under the terms of the MIT License.

"""

import logging

__author__ = "James Glasbrenner"
__copyright__ = "Copyright 2017, Mason DataMaterials Group"
__maintainer__ = "James Glasbrenner"
__email__ = "jglasbr2@gmu.edu"
__date__ = "October 6, 2017"

logger = logging.getLogger(__name__)


class SpynsStructure(object):
    """Crystal structure defining a material system."""

    def __init__(self):
        """Build SpynsStructure directly from a pymatgen Structure object.

        Parameters
        ----------
        structure : `pymatgen.core.structure.Structure` instance
            Encodes the structural information of a material via pymatgen's
            Structure class.

        """
        pass
