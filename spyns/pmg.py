# -*- coding : utf - 8 -*-
"""Spyns interface with pymatgen.

Copyright (c) Mason DataMaterials Group.
Distributed under the terms of the MIT License.

"""

import copy
import logging
import math
import operator
import sys

import numpy as np
import pandas as pd
from pymatgen import Lattice, Structure
from pymatgen.transformations.standard_transformations import \
    SupercellTransformation

from spyns.utils import convert_spherical_to_cartesian, groupby_list

__author__ = "James Glasbrenner"
__copyright__ = "Copyright 2017, Mason DataMaterials Group"
__maintainer__ = "James Glasbrenner"
__email__ = "jglasbr2@gmu.edu"
__date__ = "October 9, 2017"

logger = logging.getLogger(__name__)


def build_structure(structure_parameters):
    """Placeholder."""
    try:
        input_method = structure_parameters.get('input')
        site_properties = structure_parameters.get('site_properties')
        transformations = structure_parameters.get('transformations')

    except AttributeError:
        logger.error("Structure parameters information is not " "dict-like.")
        raise

    if input_method == "spacegroup":
        try:
            structure = from_spacegroup_to_structure(
                sg=structure_parameters["sg"],
                crystal_system=(
                    structure_parameters["lattice"]["crystal_system"]),
                lattice_parameters=(
                    structure_parameters["lattice"]["parameters"]),
                species=structure_parameters["species"],
                coords=structure_parameters["coords"],
                site_properties=site_properties,
            )

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
                format_site_properties(site_properties)
                logger.debug("Site-properties post-formatting: "
                             "%s", site_properties)
                set_site_properties(structure, site_properties)

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
        structure = apply_structure_transformations(structure, transformations)

    structure = sort_spin_sites(pmg_structure=structure)

    return structure


def from_spacegroup_to_structure(sg, crystal_system, lattice_parameters,
                                 species, coords, site_properties=None):
    """Placeholder."""
    lattice = getattr(Lattice, crystal_system)

    try:
        logger.debug("Site-properties pre-formatting: " "%s", site_properties)
        format_site_properties(site_properties)
        logger.debug("Site-properties post-formatting: " "%s", site_properties)

    except AttributeError:
        logger.debug("site_properties is not a dict-like object. "
                     "Structure object will not have site properties.")

    structure = Structure.from_spacegroup(
        sg=sg,
        lattice=lattice(**lattice_parameters),
        species=species,
        coords=coords,
        site_properties=site_properties,
    )

    return structure


def set_site_properties(structure, site_properties):
    """Placeholder."""
    for property_name, site_values in site_properties.items():
        structure.add_site_property(property_name, site_values)


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
                logging.error("Unsupported property name %s", property_name)
                sys.exit("Unsupported property name.")

        # Sort the list by site index
        sorted(list_of_site_values, key=operator.itemgetter(0))
        logger.debug("Site values for %s property: "
                     "%s", property_name, list_of_site_values)

        # Overwrite property_name dict with structured list
        site_properties[property_name] = [
            site[1] for site in list_of_site_values
        ]


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


def sort_spin_sites(pmg_structure):
    """Placeholder."""
    sorted_pmg_structure = pmg_structure.get_sorted_structure(
        lambda x: (x.frac_coords[2], x.frac_coords[1], x.frac_coords[0]))

    return sorted_pmg_structure


def find_all_neighbors(pmg_structure, cutoff):
    """Placeholder."""
    all_neighbors = pmg_structure.get_all_neighbors(r=cutoff,
                                                    include_index=True)
    grouped_neighbors = group_all_neighbors(pmg_structure, all_neighbors)

    return grouped_neighbors


def group_all_neighbors(pmg_structure, all_neighbors):
    """Placeholder."""
    sublattices_set = set(pmg_structure.site_properties["sublattice"])
    unique_sublattice_distances = {}
    for sublattice_id in sublattices_set:
        unique_sublattice_distances["{0}".format(sublattice_id)] = {
            "{0}".format(x): []
            for x in sublattices_set
        }
    logger.debug("Pre-allocated dict for sublattice distances = %s",
                 unique_sublattice_distances)

    for site_index, site_neighbors in enumerate(all_neighbors):
        site_sublattice_id = pmg_structure[site_index].properties["sublattice"]
        unique_site_distances = unique_sublattice_distances["{0}".format(
            site_sublattice_id)]

        for neighbor in site_neighbors:
            neighbor_sublattice_id = neighbor[0].properties["sublattice"]
            neighbor_distance = neighbor[1]
            unique_neighbor_distances = unique_site_distances["{0}".format(
                neighbor_sublattice_id)]
            new_distance_test = [
                math.isclose(neighbor_distance, x)
                for x in unique_neighbor_distances
            ]
            if not any(new_distance_test):
                unique_neighbor_distances.append(neighbor_distance)

    # Sort unique distances between sublattice-sublattice pairs
    for sublattice_pairs in unique_sublattice_distances.values():
        for unique_distances in sublattice_pairs.values():
            unique_distances.sort()

    all_neighbors_sorted = []
    for site_index, site_neighbors in enumerate(all_neighbors):
        site_sublattice_id = pmg_structure[site_index].properties["sublattice"]
        unique_site_distances = unique_sublattice_distances["{0}".format(
            site_sublattice_id)]

        site_neighbors_numbered = []
        for neighbor in site_neighbors:
            neighbor_sublattice_id = neighbor[0].properties["sublattice"]
            neighbor_distance = neighbor[1]
            unique_neighbor_distances = unique_site_distances["{0}".format(
                neighbor_sublattice_id)]
            for neighbor_number, test_distance in enumerate(
                    unique_neighbor_distances, start=1):
                if math.isclose(neighbor_distance, test_distance):
                    update_row = [x for x in neighbor]
                    update_row[1] = tuple([update_row[1], neighbor_number])
            site_neighbors_numbered.append(tuple(update_row))

        site_neighbors_sorted = sorted(
            site_neighbors_numbered,
            key=lambda x: (x[0].properties["sublattice"], x[1][1]))

        all_neighbors_sorted.append(site_neighbors_sorted)

    all_neighbors_reduced = []
    for site_neighbors in all_neighbors_sorted:
        grouped_neighbors = groupby_list(
            iterable=site_neighbors,
            key=lambda x: (x[0].properties["sublattice"], x[1][1]))
        site_neighbors_reduced = []
        for neighbor in grouped_neighbors:
            neighbor_sublattice_id = neighbor[0][0]
            neighbor_number = neighbor[0][1]
            neighbor_list = neighbor[1]
            neighbor_indices = tuple(x[2] for x in neighbor_list)
            site_neighbors_reduced.append((neighbor_sublattice_id,
                                           neighbor_number, neighbor_indices))
        all_neighbors_reduced.append(site_neighbors_reduced)

    build_dataframe = []
    for site_index, site_neighbors in enumerate(all_neighbors_reduced):
        site_sublattice_id = pmg_structure[site_index].properties["sublattice"]
        for grouped_neighbors in site_neighbors:
            neighbor_sublattice_id = grouped_neighbors[0]
            neighbor_number = grouped_neighbors[1]
            neighbor_indices = grouped_neighbors[2]
            for neighbor_index in neighbor_indices:
                build_dataframe.append(
                    tuple((site_index, neighbor_index, site_sublattice_id,
                           neighbor_sublattice_id, neighbor_number)))
    neighbors_df = np.array(build_dataframe, dtype={
        "names": [
            "site_i", "site_j", "sublattice_i", "sublattice_j",
            "neighbor_number"
        ],
        "formats": ["i8", "i8", "i8", "i8", "i8"]
    })
    neighbors_df = pd.DataFrame(data=neighbors_df)

    build_dataframe = []
    for (site_sublattice_id,
         neighbor_sublattice_distances) in unique_sublattice_distances.items():
        for (neighbor_sublattice_id,
             neighbor_distances) in neighbor_sublattice_distances.items():
            for neighbor_number, distance in enumerate(neighbor_distances,
                                                       start=1):
                build_dataframe.append(
                    tuple((site_sublattice_id, neighbor_sublattice_id,
                           neighbor_number, distance)))
    distances_df = np.array(build_dataframe, dtype={
        "names":
        ["sublattice_i", "sublattice_j", "neighbor_number", "distance"],
        "formats": ["i8", "i8", "i8", "f8"]
    })
    distances_df = pd.DataFrame(data=distances_df)

    return neighbors_df, distances_df
