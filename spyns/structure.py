# -*- coding: utf-8 -*-
"""Hamiltonian constructor for the Ising model.

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


class Hamiltonian(object):
    """Pair-wise interactions between Ising spins."""

    def __init__(self):
        """Build Hamiltonian object using site and neighbor information.

        Parameters
        ----------
        sites : list
            Positions of `n` sites indexed from 0 to `n` - 1

        neighbors : dict
            Each site's neighbors grouped and sorted by distance and identified
            by site index

        """
        pass
