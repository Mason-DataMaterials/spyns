# -*- coding : utf - 8 -*-
"""Data structures for Spyns simulations.

Copyright (c) Mason DataMaterials Group.
Distributed under the terms of the MIT License.

"""
import copy
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pymatgen
import ruamel.yaml
from monty.json import MSONable
from pymatgen.alchemy.materials import TransformedStructure
from pymatgen.transformations.standard_transformations import \
    SupercellTransformation
from ruamel.yaml.scanner import ScannerError

from spyns.pmg import build_structure, find_all_neighbors
from spyns.utils import convert_spherical_to_cartesian

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
        self._neighbor_distances = None
        self._spins = None
        self._sublattices = None
        self._exchange = None
        self._energy = None
        self._totenergy = None
        self._scaling_factors = {}

        if pmg_structure and isinstance(pmg_structure, pymatgen.Structure):
            logger.debug("Input structure appears to be a pymatgen.Structure "
                         "object.")
            self.pmg_structure = TransformedStructure(structure=pmg_structure)

        else:
            logger.debug("Input structure appears to not be a "
                         "pymatgen.Structure object. Checking if dict-like "
                         "object.")

            try:
                pmg_structure.keys()
                self.pmg_structure = TransformedStructure(
                    structure=self._build_structure(pmg_structure))

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

    def scale_supercell(self, scale_abc):
        """Placeholder."""
        scale_a = scale_abc[0]
        scale_b = scale_abc[1]
        scale_c = scale_abc[2]

        scale_transformation = SupercellTransformation.from_scaling_factors(
            scale_a=scale_a, scale_b=scale_b, scale_c=scale_c)
        self.pmg_structure.append_transformation(scale_transformation)

    @property
    def spins(self):
        """Placeholder."""
        return self._spins

    @spins.setter
    def spins(self, spins_input):
        """Placeholder."""
        build_dataframe = []
        if type(spins_input) is list or type(spins_input) is np.ndarray:
            for site_index, spin in enumerate(spins_input):
                spin_x = 0
                spin_y = 0
                spin_z = spin
                build_dataframe.append(
                    tuple([site_index, spin_x, spin_y, spin_z]))

        elif type(spins_input) is dict:
            for site_index, spin_components in spins_input.items():
                rho = spin_components["rho"]
                phi = spin_components["phi"]
                theta = spin_components["theta"]
                spin_x, spin_y, spin_z = convert_spherical_to_cartesian(
                    rho=rho, phi=phi, theta=theta)
                build_dataframe.append([site_index, spin_x, spin_y, spin_z])

        spins_df = np.array(build_dataframe, dtype={
            "names": ["site_i", "spin_x", "spin_y", "spin_z"],
            "formats": ["i8", "f8", "f8", "f8"]
        })
        spins_df = pd.DataFrame(data=spins_df)
        self._spins = spins_df

    @property
    def neighbors(self):
        """Placeholder."""
        return self._neighbors

    @neighbors.setter
    def neighbors(self, cutoff):
        """Placeholder."""
        pmg_structure = self.pmg_structure.final_structure
        structure_history_id = len(self.pmg_structure.structures) - 1
        neighbors_df, distances_df = find_all_neighbors(
            pmg_structure=pmg_structure, cutoff=cutoff)

        if self._neighbors:
            self._neighbors = pd.concat([
                self._neighbors,
                neighbors_df.assign(history_id=structure_history_id)
            ], ignore_index=True)

        else:
            self._neighbors = neighbors_df.assign(
                history_id=structure_history_id)

        if self._sublattices:
            self._sublattices = pd.merge(
                left=self._sublattices, right=distances_df,
                on=["sublattice_i", "sublattice_j", "neighbor_number"])

        else:
            self._sublattices = distances_df

    @property
    def pairs_df(self):
        """Placeholder."""
        return self._sublattices

    @property
    def sublattices(self):
        """Placeholder."""
        return None

    @sublattices.setter
    def sublattices(self, value):
        """Placeholder."""
        self.pmg_structure.final_structure.add_site_property(
            "sublattice", value)

    @property
    def exchange(self):
        """Placeholder."""
        return self._exchange

    @exchange.setter
    def exchange(self, exchange_list):
        """Placeholder."""
        build_dataframe = []
        for exchange_ij in exchange_list:
            sublattice_i, sublattice_j, neighbor_number, exchange = exchange_ij
            build_dataframe.append(
                tuple([sublattice_i, sublattice_j, neighbor_number, exchange]))
        exchange_df = np.array(build_dataframe, dtype={
            "names":
            ["sublattice_i", "sublattice_j", "neighbor_number", "exchange_ij"],
            "formats": ["i8", "i8", "i8", "f8"]
        })
        exchange_df = pd.DataFrame(data=exchange_df)
        if self._sublattices is not None:
            self._sublattices = pd.merge(
                left=self._sublattices, right=exchange_df,
                on=["sublattice_i", "sublattice_j", "neighbor_number"])
        else:
            self._sublattices = exchange_df


class SimulationData(MSONable):
    """Placeholder."""

    def __init__(self):
        """Placeholder."""
        pass
