# -*- coding: utf-8 -*-
"""Monte Carlo simulations of the Ising model using the Metropolis algorithm.

Copyright (c) Mason DataMaterials Group.
Distributed under the terms of the MIT License.

Notes
-----
Our code is an adaptation of the Ising model code example in chapter 7.7.1,
pages 318–319, of *Annotated Algorithms in Python* [ann_algo_python]_,
which was released under the BSDv3 License. The starting code remains under
that license, while all subsequent changes and new features are released under
the MIT license.

References
----------
.. [ann_algo_python] Massimo Di Pierro, *Annotated Algorithms in Python: With
   Applications in Physics, Biology, and Finance*, 1st ed. (Experts4solutions,
   Lexington, KY, 2014).

"""

from pathlib import Path
import logging
import math
import random

import numpy as np
import pandas as pd

__author__ = "Swabir Silayi, James Glasbrenner"
__copyright__ = "Copyright 2017, Mason DataMaterials Group"
__maintainer__ = "James Glasbrenner"
__email__ = "jglasbr2@gmu.edu"
__date__ = "August 31, 2017"

logger = logging.getLogger(__name__)


class Ising(object):
    """Monte Carlo simulation of the Ising model."""

    coords = []

    def __init__(self, number_sites_along_xyz):
        """Set the number of sites in the cubic lattice.

        Parameters
        ----------
        nsites_along_xyz : int
            Sets the number of sites along each side of the cubic lattice.
            Defines a lattice with total sites equal to nsites_along_xyz**3

        """
        self.number_sites_along_xyz = number_sites_along_xyz
        self.site_spin = [
            [
                [1 for x in range(number_sites_along_xyz)]
                for y in range(number_sites_along_xyz)
            ] for z in range(number_sites_along_xyz)
        ]
        self.magnetization = number_sites_along_xyz**3

    def __getitem__(self, site_index):
        """Return the spin value at the specified site index.

        Parameters
        ----------
        site_index : int
            A site index from 0 to number_sites_along_xyz**3

        """
        number_sites_along_xyz = self.number_sites_along_xyz
        site_x, site_y, site_z = site_index

        index_x = (site_x + number_sites_along_xyz) % number_sites_along_xyz
        index_y = (site_y + number_sites_along_xyz) % number_sites_along_xyz
        index_z = (site_z + number_sites_along_xyz) % number_sites_along_xyz

        return self.site_spin[index_x][index_y][index_z]

    def __setitem__(self, site_index, new_spin_value):
        """Set the spin value at the specified site index.

        Parameters
        ----------
        site_index : int
            A site index from 0 to number_sites_along_xyz**3

        new_spin_value : int
            In the Ising model, the spin value is restricted to values of 1
            or -1

        """
        number_sites_along_xyz = self.number_sites_along_xyz
        site_x, site_y, site_z = site_index

        index_x = (site_x + number_sites_along_xyz) % number_sites_along_xyz
        index_y = (site_y + number_sites_along_xyz) % number_sites_along_xyz
        index_z = (site_z + number_sites_along_xyz) % number_sites_along_xyz

        self.site_spin[index_x][index_y][index_z] = new_spin_value

    def run_simulation(self, steps, temperature, external_field):
        """Run Monte Carlo simulation.

        Parameters
        ----------
        steps : int
            Number of time steps to run the simulation.

        external_field : float
            External magnetic field in units of (check units)

        temperature : float
            System temperature in units of (check units)

        Returns
        -------
        simulation_statistics
            The mean and standard deviation of the magnetization at a given
            temperature.

        """
        magnetization_history = []

        for _ in range(steps):
            self.step(
                temperature=float(temperature),
                external_field=float(external_field),
            )
            magnetization_history.append(self.magnetization)

        logger.info(msg="h={0}, t={1}".format(external_field, temperature))

        mean_magnetization = np.mean(magnetization_history)
        std_dev_magnetization = np.std(magnetization_history)
        simulation_statistics = (
            temperature,
            mean_magnetization,
            std_dev_magnetization,
        )

        self.write_state_snapshot_to_disk(
            temperature=temperature,
            external_field=external_field,
        )

        return simulation_statistics

    def step(self, external_field, temperature):
        """Advance the Monte Carlo simulation by one time step.

        Parameters
        ----------
        external_field : float
            External magnetic field in units of (check units)

        temperature : float
            System temperature in units of (check units)

        """
        site_x, site_y, site_z = self._select_spin_site()
        neighbors = self._find_site_neighbors(site_x, site_y, site_z)
        change_in_total_energy = self._calculate_spin_flip_energy(
            site_x,
            site_y,
            site_z,
            neighbors,
            external_field,
        )
        self._keep_or_reject_spin_flip(
            change_in_total_energy=change_in_total_energy,
            temperature=temperature,
            site_x=site_x,
            site_y=site_y,
            site_z=site_z,
        )

    def write_state_snapshot_to_disk(self, external_field, temperature):
        """Save system state to output file.

        Parameters
        ----------
        external_field : float
            External magnetic field in units of (check units)

        temperature : float
            System temperature in units of (check units)

        """
        output_directory = (
            Path(".").cwd() /
            "results/h_{0}/configurations".format(external_field)
        )
        output_directory.mkdir(parents=True, exist_ok=True)
        output_filename = output_directory / "t_{0}.csv".format(temperature)

        state_full_system = []
        column_names = ["site_x", "site_y", "site_z", "spin"]
        self._gather_current_state(state_full_system)

        _save_dataset_to_disk(
            dataset=state_full_system,
            column_names=column_names,
            output_filename=output_filename,
        )

    def _calculate_spin_flip_energy(
            self,
            site_x,
            site_y,
            site_z,
            neighbors,
            external_field,
    ):

        site_spin = self[site_x, site_y, site_z]
        neighbor_spins = []
        for neighbor_x, neighbor_y, neighbor_z in neighbors:
            neighbor_spins.append(self[neighbor_x, neighbor_y, neighbor_z])

        change_in_total_energy = (
            -2.0 * site_spin * (external_field + sum(neighbor_spins))
        )

        return change_in_total_energy

    def _gather_current_state(self, state_full_system):

        for i in range(self.number_sites_along_xyz):
            for j in range(self.number_sites_along_xyz):
                for k in range(self.number_sites_along_xyz):
                    state_single_site = np.array(
                        [i, j, k, self.site_spin[i][j][k]]
                    )
                    state_full_system.append(state_single_site)

    def _keep_or_reject_spin_flip(
            self,
            change_in_total_energy,
            temperature,
            site_x,
            site_y,
            site_z,
    ):

        if change_in_total_energy > temperature * math.log(random.random()):

            self[site_x, site_y, site_z] = -self[site_x, site_y, site_z]
            self.magnetization += 2 * self[site_x, site_y, site_z]

    def _select_spin_site(self):

        site_x, site_y, site_z = (
            random.randint(0, self.number_sites_along_xyz - 1),
            random.randint(0, self.number_sites_along_xyz - 1),
            random.randint(0, self.number_sites_along_xyz - 1),
        )

        return site_x, site_y, site_z

    @staticmethod
    def _find_site_neighbors(site_x, site_y, site_z):

        neighbors = [
            (site_x - 1, site_y, site_z),
            (site_x + 1, site_y, site_z),
            (site_x, site_y - 1, site_z),
            (site_x, site_y - 1, site_z),
            (site_x, site_y, site_z - 1),
            (site_x, site_y, site_z + 1),
        ]

        return neighbors


def main(
        number_sites_along_xyz=10,
        steps=25000,
        external_field_sweep_start=1,
        external_field_sweep_end=11,
        temperature_sweep_start=1,
        temperature_sweep_end=11,
):
    """Run simulation over a sweep of temperature and external field values.

    Parameters
    ----------
    number_sites_along_xyz : int
    steps : int
    external_field_sweep_start : int
    external_field_sweep_end : int
    temperature_sweep_start : int
    temperature_sweep_end : int

    """
    ising = Ising(number_sites_along_xyz=number_sites_along_xyz)
    column_names = [
        "temperature",
        "mean_magnetization",
        "std_dev_magnetization",
    ]

    _external_field_sweep(
        ising=ising,
        steps=steps,
        column_names=column_names,
        external_field_sweep_start=external_field_sweep_start,
        external_field_sweep_end=external_field_sweep_end,
        temperature_sweep_start=temperature_sweep_start,
        temperature_sweep_end=temperature_sweep_end,
    )


def _external_field_sweep(
        ising,
        steps,
        column_names,
        external_field_sweep_start,
        external_field_sweep_end,
        temperature_sweep_start,
        temperature_sweep_end,
):
    for external_field in range(external_field_sweep_start,
                                external_field_sweep_end):

        simulation_results = []

        output_directory = (
            Path(".").cwd() / "results/h_{0}".format(external_field)
        )
        output_directory.mkdir(parents=True, exist_ok=True)
        output_filename = (output_directory / "mg.csv")

        _temperature_sweep(
            ising=ising,
            steps=steps,
            external_field=external_field,
            simulation_results=simulation_results,
            sweep_start=temperature_sweep_start,
            sweep_end=temperature_sweep_end,
        )

        _save_dataset_to_disk(
            dataset=simulation_results,
            column_names=column_names,
            output_filename=output_filename,
        )


def _save_dataset_to_disk(dataset, column_names, output_filename):

    db = pd.DataFrame(
        data=np.array(dataset),
        columns=column_names,
    )

    db.to_csv(
        path_or_buf="{0}".format(output_filename),
        index=False,
    )


def _temperature_sweep(
        ising,
        steps,
        external_field,
        simulation_results,
        sweep_start,
        sweep_end,
):
    for temperature in range(sweep_start, sweep_end):

        simulation_statistics = ising.run_simulation(
            steps=steps,
            temperature=temperature,
            external_field=external_field,
        )
        simulation_results.append(simulation_statistics)


if __name__ == '__main__':

    number_sites_along_xyz = 10
    steps = 25000

    main(
        number_sites_along_xyz=number_sites_along_xyz,
        steps=steps,
    )
