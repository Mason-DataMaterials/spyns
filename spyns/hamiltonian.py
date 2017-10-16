# -*- coding: utf-8 -*-
"""Hamiltonian constructor for the pairwise spin interaction model.

Copyright (c) Mason DataMaterials Group.
Distributed under the terms of the MIT License.

"""

import logging
import sys

import numpy as np
import pandas as pd

from spyns.utils import convert_spherical_to_cartesian

__author__ = "James Glasbrenner"
__copyright__ = "Copyright 2017, Mason DataMaterials Group"
__maintainer__ = "James Glasbrenner"
__email__ = "jglasbr2@gmu.edu"
__date__ = "October 6, 2017"

logger = logging.getLogger(__name__)


class Hamiltonian(object):
    """Pair-wise interactions between Ising spins."""

    def __init__(self, pairwise_neighbors_df, pairwise_interactions,
                 pairwise_interaction_type="ising"):
        """Build Hamiltonian object using site and neighbor information.

        Parameters
        ----------
        sites : list
            Positions of `n` sites indexed from 0 to `n` - 1

        neighbors : dict
            Each site's neighbors grouped and sorted by distance and identified
            by site index

        """
        if (pairwise_interaction_type == "ising"
                or pairwise_interaction_type == "heisenberg"):
            self.pairwise_interaction_type = pairwise_interaction_type

        else:
            logger.error("Unrecognized interaction type: %s",
                         pairwise_interaction_type)
            sys.exit(1)

        self._total_energy = 0

        build_dataframe = []
        for (interaction_index,
             interaction_parameter) in pairwise_interactions.items():
            sublattice_i, sublattice_j, neighbor_number = interaction_index
            interaction_coefficient = interaction_parameter
            if sublattice_i == sublattice_j:
                build_dataframe.append(
                    tuple([
                        sublattice_i, sublattice_j, neighbor_number,
                        interaction_coefficient
                    ]))
            else:
                build_dataframe.append(
                    tuple([
                        sublattice_i, sublattice_j, neighbor_number,
                        interaction_coefficient
                    ]))
                build_dataframe.append(
                    tuple([
                        sublattice_j, sublattice_i, neighbor_number,
                        interaction_coefficient
                    ]))
        pairwise_interactions_df = np.array(build_dataframe, dtype={
            "names": [
                "sublattice_i", "sublattice_j", "neighbor_number",
                "interaction_ij"
            ],
            "formats": ["i8", "i8", "i8", "f8"]
        })
        pairwise_interactions_df = pd.DataFrame(data=pairwise_interactions_df)

        self.pairwise_interactions_df = pairwise_interactions_df
        self.pairwise_neighbors_df = pairwise_neighbors_df.sort_values(
            by=[
                "site_i", "sublattice_i", "sublattice_j", "neighbor_number",
                "site_j"
            ])

        pairwise_neighbors_df = (self.pairwise_neighbors_df.merge(
            right=self.pairwise_interactions_df,
            on=["sublattice_i", "sublattice_j",
                "neighbor_number"]).sort_values(by=[
                    "site_i", "sublattice_i", "sublattice_j",
                    "neighbor_number", "site_j"
                ]))
        self._exchange_array = (
            pairwise_neighbors_df["interaction_ij"].as_matrix())
        self._sites_array = pairwise_neighbors_df["site_i"].as_matrix()
        self._neighbors_array = pairwise_neighbors_df["site_j"].as_matrix()

        neighbor_lookup_df = (self.pairwise_neighbors_df.groupby(
            by=["site_i", "sublattice_i", "sublattice_j", "neighbor_number"],
            as_index=False).agg({
                "site_j": "count"
            }).rename(columns={"site_j": "neighbor_count"})).sort_values([
                "site_i", "sublattice_i", "sublattice_j", "neighbor_number"
            ]).reset_index(drop=True).groupby("site_i", as_index=False).agg({
                "neighbor_count":
                "sum"
            })
        self._neighbors_cumindex_array = (
            neighbor_lookup_df["neighbor_count"].cumsum().as_matrix())

    @property
    def exchange(self):
        """Placeholder."""
        return self._exchange_array

    @property
    def neighbors(self):
        """Placeholder."""
        return self._neighbors_array

    @property
    def neighbors_cumindex(self):
        """Placeholder."""
        return self._neighbors_cumindex_array

    @property
    def spins(self):
        """Placeholder."""
        self.spins_df.sort_values(by="site_i", inplace=True)

        if self.pairwise_interaction_type == "ising":
            spins_array = self.spins_df["spin_z"].as_matrix()

        elif self.pairwise_interaction_type == "heisenberg":
            spins_array = self.spins_df[["spin_x", "spin_y",
                                         "spin_z"]].as_matrix().flatten()
        else:
            logger.error("Unrecognized interaction type: %s",
                         self.pairwise_interaction_type)
            sys.exit(1)

        return spins_array

    @spins.setter
    def spins(self, spins_state):
        """Placeholder."""
        build_dataframe = []
        if type(spins_state) is list or type(spins_state) is np.ndarray:
            for site_index, spin in enumerate(spins_state):
                spin_x = 0
                spin_y = 0
                spin_z = spin
                build_dataframe.append(
                    tuple([site_index, spin_x, spin_y, spin_z]))

        elif type(spins_state) is dict:
            for site_index, spin_components in spins_state.items():
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
        self.spins_df = spins_df

        nsites = len(self.spins_df)
        if self.pairwise_interaction_type == "ising":
            spins_array = self.spins_df["spin_z"].as_matrix()
            self._total_energy = 0.5 * (
                self._exchange_array * spins_array[self._sites_array] *
                spins_array[self._neighbors_array.tolist()]).sum() / nsites

        elif self.pairwise_interaction_type == "heisenberg":
            spins_array = self.spins_df[["spin_x", "spin_y",
                                         "spin_z"]].as_matrix().flatten(
                                             order="F")
            heisenberg_sites_array = []
            for x in range(3):
                tmp_array = self._sites_array + x * nsites
                heisenberg_sites_array.extend(tmp_array.tolist())
            heisenberg_neighbors_array = []
            for x in range(3):
                tmp_array = self._neighbors_array + x * nsites
                heisenberg_neighbors_array.extend(tmp_array.tolist())
            heisenberg_exchange_array = np.repeat(self._exchange_array, 3)

            self._total_energy = (0.5 * (
                heisenberg_exchange_array * spins_array[heisenberg_sites_array]
                * spins_array[heisenberg_neighbors_array]).sum() / nsites)

        else:
            logger.error("Unrecognized interaction type: %s",
                         self.pairwise_interaction_type)
            sys.exit(1)

    @property
    def energy(self):
        """Placeholder."""
        return self._total_energy
