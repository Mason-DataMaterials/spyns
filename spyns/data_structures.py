# -*- coding : utf - 8 -*-
"""Data structures for Spyns simulations.

Copyright (c) Mason DataMaterials Group.
Distributed under the terms of the MIT License.

"""

import copy
import logging
from pathlib import Path

import pymatgen
import ruamel.yaml
from monty.json import MSONable
from ruamel.yaml.scanner import ScannerError

from spyns.pmg import build_structure

__author__ = "James Glasbrenner"
__copyright__ = "Copyright 2017, Mason DataMaterials Group"
__maintainer__ = "James Glasbrenner"
__email__ = "jglasbr2@gmu.edu"
__date__ = "October 6, 2017"

logger = logging.getLogger(__name__)


class SpynsSystem(MSONable):
    """Crystal structure and parameters defining a material system."""

    def __init__(self, pmg_structure):
        """Build SpynsSystem directly from a pymatgen Structure object.

        Parameters
        ----------
        structure : `pymatgen.core.structure.Structure` instance
            Encodes the structural information of a material via pymatgen's
            Structure class.

        """
        self._neighbors = None
        self._spins = None
        self._sublattices = None

        if pmg_structure and isinstance(pmg_structure, pymatgen.Structure):
            logger.debug("Input structure appears to be a pymatgen.Structure "
                         "object.")
            self.pmg_structure = pmg_structure

        else:
            logger.debug("Input structure appears to not be a "
                         "pymatgen.Structure object. Checking if dict-like "
                         "object.")

            try:
                pmg_structure.keys()
                self.pmg_structure = self._build_structure(pmg_structure)

            except AttributeError:
                logger.debug("Input structure is not a dict-like object. "
                             "It is instead the type %s. ",
                             type(pmg_structure))
                raise

            except TypeError:
                logger.debug("Input structure is not a dict-like object. "
                             "It is instead the type %s. ",
                             type(pmg_structure))
                raise

            except NameError:
                logger.debug("Input structure is not defined.")
                raise

    @classmethod
    def from_yaml(cls, filename: str):
        """Create SpynsSystem by reading parameters from a YAML file.

        The YAML file defines a crystal structure by pointing to a VASP POSCAR
        file or CIF file, or by providing spacegroup information. If using
        a POSCAR or CIF file, schematically the YAML file should contain:

        .. code-block:: yaml
           input: file
           filepath: relative or absolute path to CIF or POSCAR file

        If using spacegroup information, then the schematic for the YAML file
        is:

        .. code-block:: yaml
           input: spacegroup
           sg: number or Hermann-Mauguin (international) notation
           lattice:
             crystal_system: cubic, tetragonal, etc.
             constants: (different for each crystal system)
               a: float
               c: float
           species:
             - species1
             - species2
           coords:
             - species1_coordinates
             - species2_coordinates

        Optional site properties allow users to include additional information
        and parameters, such as sublattice labels or what to use for the
        initial spin configuration. The following schematic for the YAML file
        defines the spin configuration and sublattice labels for a two-species
        system:

        .. code-block:: yaml
           site_properties:
             magmom:
               1:
                 rho: species1_magmom_rho
                 phi: species1_magmom_phi
                 theta: species1_magmom_theta
               2:
                 x: species2_magmom_x
                 y: species2_magmom_y
                 z: species2_magmom_z
             sublattice:
               1: species1_sublattice_numerical_id
               2: species2_sublattice_numerical_id

        If you are specifying site properties for a structure generated using
        a POSCAR or CIF file, then for each site property you need to provide
        a value for each site, from 1 to N, in the file. If instead you
        provide spacegroup information, then the range 1 to N indexes the
        species list.

        Supercell transformations are also optional, but can be specified to
        build out structures. This becomes useful when you need to define
        sublattices. The schematic for defining supercell transformations is:

        .. code-block:: yaml
           transformations:
             primitive: bool
             matrix:
               - - 1
                 - -1
                 - 0
               - - 1
                 - 1
                 - 0
               - - 0
                 - 0
                 - 1
             scale:
               a: 1
               b: 1
               c: 1

        Setting the ``primitive`` option to `True` transforms the structure
        to its primitive cell, which is applied first. The ``matrix`` option
        defines the scaling matrix, which allows you to rotate and scale
        the vectors defining the supercell, and is applied second. The `scale`
        option is a (simpler) alternative to the matrix option and is applied
        last, allowing you to extend the supercell along the ``a``, ``b``,
        and/or ``c`` supercell parameters. You can specify one, two, or all
        three options in your file, depending on your preferences.

        Parameters
        ----------
        filename : str
            Path to YAML file defining the system's crystal structure and
            chemical composition.

        Returns
        -------
        cls
            An instance of the SpynsSystem class

        """
        try:
            logger.debug("Reading file = %s", filename)
            yaml = ruamel.yaml.YAML()
            yaml_dict = yaml.load(Path(filename))
            spyns_structure = build_structure(yaml_dict)
            logger.debug("spyns_structure = %s", spyns_structure)

        except TypeError:
            logger.debug("Reading dict-like object = %s", filename)
            yaml_dict = copy.deepcopy(filename)
            spyns_structure = build_structure(yaml_dict)
            logger.debug("spyns_structure = %s", spyns_structure)

        except ScannerError:
            logger.error("Input is not a YAML file or dict-like object!")
            raise

        return cls(pmg_structure=spyns_structure)

    @property
    def spins(self):
        """Placeholder."""
        return self._spins

    @spins.setter
    def spins(self, value):
        """Placeholder."""
        pass

    @property
    def neighbors(self):
        """Placeholder."""
        return self._neighbors

    @neighbors.setter
    def neighbors(self, cutoff):
        """Placeholder."""
        pmg_neighbors = self.pmg_structure.get_all_neighbors(r=cutoff)
        self._neighbors = pmg_neighbors

    @property
    def sublattices(self):
        """Placeholder."""
        return self._sublattices

    @sublattices.setter
    def sublattices(self, value):
        """Placeholder."""
        pass


class SimulationData(MSONable):
    """Placeholder."""

    def __init__(self):
        """Placeholder."""
        pass
