# -*- coding : utf - 8 -*-
"""Data structures for Spyns simulations.

Copyright (c) Mason DataMaterials Group.
Distributed under the terms of the MIT License.

"""

import copy
import logging
import operator
import sys
from pathlib import Path

import ruamel.yaml
from monty.json import MSONable
from pymatgen import Lattice, Structure
from pymatgen.transformations.standard_transformations import \
    SupercellTransformation
from ruamel.yaml.scanner import ScannerError

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
        self._spins = None
        self._sublattices = None

        if pmg_structure and isinstance(pmg_structure, Structure):
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
            spyns_structure = cls.build_structure(yaml_dict)
            logger.debug("spyns_structure = %s", spyns_structure)

        except TypeError:
            logger.debug("Reading dict-like object = %s", filename)
            yaml_dict = copy.deepcopy(filename)
            spyns_structure = cls.build_structure(yaml_dict)
            logger.debug("spyns_structure = %s", spyns_structure)

        except ScannerError:
            logger.error("Input is not a YAML file or dict-like object!")
            raise

        return cls(pmg_structure=spyns_structure)

    @staticmethod
    def build_structure(structure_parameters):
        """Placeholder."""
        cls = SpynsSystem

        try:
            input_method = structure_parameters.get('input')
            site_properties = structure_parameters.get('site_properties')
            transformations = structure_parameters.get('transformations')

        except AttributeError:
            logger.error("Structure parameters information is not "
                         "dict-like.")
            raise

        if input_method == "spacegroup":
            try:
                structure = cls.from_spacegroup_to_structure(
                    sg=structure_parameters["sg"],
                    crystal_system=(
                        structure_parameters["lattice"]["crystal_system"]),
                    lattice_parameters=(
                        structure_parameters["lattice"]["parameters"]),
                    species=structure_parameters["species"],
                    coords=structure_parameters["coords"],
                    site_properties=site_properties, )

            except TypeError:
                logger.error("At least one of the structure spacegroup "
                             "parameters is an incorrect type.", exc_info=True)
                raise

            except KeyError:
                logger.error("At least one of the structure spacegroup "
                             "parameters is missing.", exc_info=True)
                raise

        elif input_method == "file":
            try:
                structure_filepath = structure_parameters["file"]
                structure = Structure.from_file(structure_filepath)

            except FileNotFoundError:
                logger.error("Input file not found.", exc_info=True)
                raise

            except KeyError:
                logger.error("The structure filepath parameter is not found.",
                             exc_info=True)
                raise

            if site_properties:
                try:
                    logger.debug("Site-properties pre-formatting: "
                                 "%s", site_properties)
                    cls.format_site_properties(site_properties)
                    logger.debug("Site-properties post-formatting: "
                                 "%s", site_properties)
                    cls.set_site_properties(structure, site_properties)

                except TypeError:
                    logger.error("At least one of the site properties is not "
                                 "structured correctly.", exc_info=True)
                    raise

                except AttributeError:
                    logger.error("Site properties need to be in a dict-like "
                                 "format.", exc_info=True)
                    raise

        else:
            logger.debug("The specified input method \"%s\" is invalid.",
                         input_method)
            sys.exit("Invalid input method")

        if transformations:
            structure = cls.apply_structure_transformations(
                structure, transformations)

        return structure

    @staticmethod
    def from_spacegroup_to_structure(sg, crystal_system, lattice_parameters,
                                     species, coords, site_properties=None):
        """Placeholder."""
        cls = SpynsSystem

        lattice = getattr(Lattice, crystal_system)

        try:
            logger.debug("Site-properties pre-formatting: "
                         "%s", site_properties)
            cls.format_site_properties(site_properties)
            logger.debug("Site-properties post-formatting: "
                         "%s", site_properties)

        except AttributeError:
            logger.debug("site_properties is not a dict-like object. "
                         "Structure object will not have site properties.")

        structure = Structure.from_spacegroup(
            sg=sg,
            lattice=lattice(**lattice_parameters),
            species=species,
            coords=coords,
            site_properties=site_properties, )

        return structure

    @staticmethod
    def set_site_properties(structure, site_properties):
        """Placeholder."""
        for property_name, site_values in site_properties.items():
            structure.add_site_property(property_name, site_values)

    @staticmethod
    def format_site_properties(site_properties):
        """Placeholder."""
        copy_of_site_properties = copy.deepcopy(site_properties)

        for property_name, sites in copy_of_site_properties.items():
            list_of_site_values = []

            for site_index, site_value in sites.items():
                if property_name == "magmom":
                    if all([
                            component_name in dict(site_value)
                            for component_name in ["x", "y", "z"]
                    ]):
                        x = site_value["x"]
                        y = site_value["y"]
                        z = site_value["z"]

                        components = [x, y, z]

                        list_of_site_values.append([site_index, components])

                    elif all([
                            component_name in dict(site_value)
                            for component_name in ["rho", "phi", "theta"]
                    ]):
                        rho = site_value["rho"]
                        phi = site_value["phi"]
                        theta = site_value["theta"]

                        components = convert_spherical_to_cartesian(
                            rho=rho, phi=phi, theta=theta)

                        list_of_site_values.append([site_index, components])

                    else:
                        logger.error("Unrecognized components for magmom.")
                        sys.exit("Unrecognized components for magmom.")

                elif property_name == "sublattice":
                    list_of_site_values.append([site_index, site_value])

                else:
                    logging.error("Unsupported property name %s",
                                  property_name)
                    sys.exit("Unsupported property name.")

            # Sort the list by site index
            sorted(list_of_site_values, key=operator.itemgetter(0))
            logger.debug("Site values for %s property: "
                         "%s", property_name, list_of_site_values)

            # Overwrite property_name dict with structured list
            site_properties[property_name] = [
                site[1] for site in list_of_site_values
            ]

    @staticmethod
    def apply_structure_transformations(structure, transformations):
        """Placeholder."""
        working_structure = copy.deepcopy(structure)
        primitive = transformations.get("primitive")
        supercell = transformations.get("supercell")

        if primitive:
            working_structure = working_structure.get_primitive_structure()

        if supercell:
            supercell_transform = supercell.get("transform")
            supercell_scaling = supercell.get("scaling")

            if supercell_transform:
                scaling_matrix = supercell_transform.get("scaling_matrix")

                transformation = (SupercellTransformation(
                    scaling_matrix=scaling_matrix))
                logger.debug("SupercellTransformation scaling matrix is "
                             "%s", transformation.scaling_matrix)
                working_structure = (
                    transformation.apply_transformation(working_structure))

            if supercell_scaling:
                scale_a = supercell_scaling["a"]
                scale_b = supercell_scaling["b"]
                scale_c = supercell_scaling["c"]

                supercell_scaling_transformation = (
                    SupercellTransformation.from_scaling_factors(
                        scale_a=scale_a, scale_b=scale_b, scale_c=scale_c))
                working_structure = (supercell_scaling_transformation.
                                     apply_transformation(working_structure))

        return working_structure

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
