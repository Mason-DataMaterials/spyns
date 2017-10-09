# -*- coding : utf - 8 -*-
"""Unit tests for loading pymatgen structures in Spyns.

Copyright (c) Mason DataMaterials Group.
Distributed under the terms of the MIT License.

"""

import logging
import logging.config
from pathlib import Path

import pytest
import ruamel.yaml

from spyns.data_structures import SpynsSystem

__author__ = "James Glasbrenner"
__copyright__ = "Copyright 2017, Mason DataMaterials Group"
__maintainer__ = "James Glasbrenner"
__email__ = "jglasbr2@gmu.edu"
__date__ = "October 9, 2017"

# Initialize yaml handler
yaml = ruamel.yaml.YAML()
yaml.default_flow_style = False

# Initialize root logger
logger_dict = yaml.load(Path("{0}".format(__file__)).parent / "logger.yaml")
logging.config.dictConfig(logger_dict)

logger = logging.getLogger(__name__)


@pytest.fixture
def pmg_iron_structure():
    """YAML-like dict that defines the crystal structure of BCC iron."""
    pmg_structure = {
        "input": "spacegroup",
        "sg": 229,
        "lattice": {
            "crystal_system": "cubic",
            "parameters": {
                "a": 2.8665,
            }
        },
        "species": ["Fe"],
        "coords": [[0.00, 0.00, 0.00]],
        "site_properties": {
            "magmom": {
                "1": {
                    "rho": 1,
                    "theta": 0,
                    "phi": 0
                }
            }
        }
    }

    return pmg_structure


def test_pmg_structure_load(pmg_iron_structure):
    """Test loading for pymatgen structure."""

    logger.debug("BEGIN UNITTEST: Load a pymatgen structure with a YAML-like "
                 "specification")
    spyns_system = SpynsSystem.from_yaml(pmg_iron_structure)
    pmg_structure = spyns_system.pmg_structure
    pmg_structure.add_site_property("sublattice", [1, 2])
    logger.debug("\nspyns_system =\n%s", pmg_structure)
    logger.debug("END UNITTEST")
