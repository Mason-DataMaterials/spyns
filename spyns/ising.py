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

from collections import Counter
from itertools import combinations_with_replacement
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

        Uses periodic boundaries for neighbors outside of grid.

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

        Uses periodic boundaries for neighbors outside of grid.

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

    def run_simulation(
            self,
            steps,
            temperature,
            external_field,
            algorithm="metropolis",
            seed=None,
    ):
        """Run Monte Carlo simulation.

        Parameters
        ----------
        steps : int
            Number of time steps to run the simulation.

        external_field : float
            External magnetic field in units of (check units)

        temperature : float
            System temperature in units of (check units)

        algorithm : str
            Sets algorithm to use for Monte Carlo simulation. Currently
            available:

            * "metropolis"
            * "slow_metropolis"

        seed : int, array_like
            Sets seed for random number generator

        Returns
        -------
        simulation_statistics
            The mean and standard deviation of the magnetization at a given
            temperature.

        """
        logger.info(
            "Starting simulation run for t={0}, h={1}, steps={2}".
            format(temperature, external_field, steps)
        )

        magnetization_history = self._monte_carlo_simulation(
            steps=steps,
            temperature=temperature,
            external_field=external_field,
            algorithm=algorithm,
            seed=seed,
        )

        mean_magnetization = np.mean(magnetization_history)
        std_dev_magnetization = np.std(magnetization_history)
        simulation_statistics = (
            temperature,
            mean_magnetization,
            std_dev_magnetization,
        )

        self._write_state_snapshot_to_disk(
            temperature=temperature,
            external_field=external_field,
        )

        return simulation_statistics

    def _gather_current_state(self, state_full_system):

        for i in range(self.number_sites_along_xyz):
            for j in range(self.number_sites_along_xyz):
                for k in range(self.number_sites_along_xyz):
                    state_single_site = np.array(
                        [i, j, k, self.site_spin[i][j][k]]
                    )
                    state_full_system.append(state_single_site)

    def _metropolis_algorithm(
            self,
            steps,
            temperature,
            external_field,
            magnetization_history,
            random_number_generator,
    ):

        random_number_array, probability_of_spin_flip = (
            self._generate_random_number_sequence(
                steps=steps,
                random_number_generator=random_number_generator,
            )
        )

        self._generate_lookup_tables(
            external_field=external_field,
            temperature=temperature,
            number_of_neighbors=6,
        )

        for step_number in range(steps):
            self._sweep_metropolis(
                temperature=temperature,
                external_field=external_field,
                random_number_array=random_number_array,
                probability_of_spin_flip=probability_of_spin_flip,
                step_number=step_number,
            )
            magnetization_history.append(self.magnetization)

    def _generate_random_number_sequence(self, steps, random_number_generator):

        random_number_array = random_number_generator.randint(
            low=0,
            high=self.number_sites_along_xyz - 1,
            size=(steps, 3),
        )
        probability_of_spin_flip = random_number_generator.uniform(
            low=0,
            high=1,
            size=steps,
        )

        return random_number_array, probability_of_spin_flip

    def _metropolis_algorithm_slow(
            self,
            steps,
            temperature,
            external_field,
            magnetization_history,
    ):

        for _ in range(steps):
            self._sweep_metropolis_slow(
                temperature=float(temperature),
                external_field=float(external_field),
            )
            magnetization_history.append(self.magnetization)

    def _monte_carlo_simulation(
            self,
            steps,
            temperature,
            external_field,
            algorithm,
            seed,
    ):
        """Run the main loop of the Markov-chain Monte Carlo simulation.

        Parameters
        ----------
        steps : int
            Number of time steps to run the simulation.

        temperature : float
            System temperature in units of (check units)

        external_field : float
            External magnetic field in units of (check units)

        algorithm : str
            Sets algorithm to use for Monte Carlo simulation. Currently
            available:

            * "metropolis"
            * "slow_metropolis"

        seed : int, array_like
            Sets seed for random number generator

        Returns
        -------
        magnetization_history : list of floats
            List of total system magnetization at each simulation time step

        """
        random_number_generator = np.random.RandomState(seed=seed)

        magnetization_history = []

        if algorithm == "metropolis":
            logger.info("Algorithm: Metropolis")
            self._metropolis_algorithm(
                steps=steps,
                temperature=temperature,
                external_field=external_field,
                magnetization_history=magnetization_history,
                random_number_generator=random_number_generator,
            )

        elif algorithm == "slow_metropolis":
            logger.info("Algorithm: Metropolis (slow implementation)")
            self._metropolis_algorithm_slow(
                steps=steps,
                temperature=temperature,
                external_field=external_field,
                magnetization_history=magnetization_history,
            )

        else:
            logger.info(
                "{0} is not a recognized algorithm. Running Metropolis "
                "algorithm by default.".format(algorithm)
            )
            self._metropolis_algorithm(
                steps=steps,
                temperature=temperature,
                external_field=external_field,
                magnetization_history=magnetization_history,
            )

        return magnetization_history

    def _sweep_metropolis(
            self,
            external_field,
            temperature,
            random_number_array,
            probability_of_spin_flip,
            step_number,
    ):
        """Single pass of Metropolis algorithm.

        Parameters
        ----------
        external_field : float
            External magnetic field in units of (check units)

        temperature : float
            System temperature in units of (check units)

        random_number_generator : np.random.RandomState object
            Container for the Mersenne Twister pseudo-random number generator

        """
        site_x, site_y, site_z = random_number_array[step_number]
        neighbors = [
            (site_x - 1, site_y, site_z),
            (site_x + 1, site_y, site_z),
            (site_x, site_y - 1, site_z),
            (site_x, site_y - 1, site_z),
            (site_x, site_y, site_z - 1),
            (site_x, site_y, site_z + 1),
        ]

        site_spin = self[site_x, site_y, site_z]
        neighbor_spins = Counter()
        for neighbor_x, neighbor_y, neighbor_z in neighbors:
            neighbor_spins["{0}".
                           format(self[neighbor_x,
                                       neighbor_y,
                                       neighbor_z,
                                       ])] += 1

        change_in_total_energy = self.delta_energy_lookup_table[
            (site_spin, neighbor_spins["1"], neighbor_spins["-1"])
        ]
        # change_in_total_energy = (
        #     -2.0 * site_spin * (external_field + neighbor_spins.sum())
        # )

        if change_in_total_energy > 0:
            self[site_x, site_y, site_z] = -self[site_x, site_y, site_z]
            self.magnetization += 2 * self[site_x, site_y, site_z]

        elif probability_of_spin_flip[step_number
                                      ] < self._find_boltzmann_probability(
                                          site_spin, neighbor_spins["1"],
                                          neighbor_spins["-1"]):
            self[site_x, site_y, site_z] = -self[site_x, site_y, site_z]
            self.magnetization += 2 * self[site_x, site_y, site_z]

    def _generate_lookup_tables(
            self,
            external_field,
            temperature,
            number_of_neighbors,
    ):

        self.delta_energy_lookup_table = {}
        self.boltzmann_lookup_table = {}
        unique_neighbor_states = combinations_with_replacement(
            (1, -1), number_of_neighbors
        )

        for neighbor_state in unique_neighbor_states:
            neighbor_spins = Counter()
            for neighbor in neighbor_state:
                neighbor_spins["{0}".format(neighbor)] += 1

            for site_spin in [1, -1]:
                self.delta_energy_lookup_table[
                    (site_spin, neighbor_spins["1"], neighbor_spins["-1"])
                ] = (
                    2.0 * site_spin * (
                        external_field +
                        (neighbor_spins["1"] - neighbor_spins["-1"])
                    )
                )

        partition_normalization = 0
        for (spin_states,
             delta_energy) in self.delta_energy_lookup_table.items():
            site_spin = spin_states[0]
            neighbors_up = spin_states[1]
            neighbors_down = spin_states[2]
            boltzmann_probability = np.exp(-delta_energy / temperature)
            self.boltzmann_lookup_table[
                (site_spin, neighbors_up, neighbors_down)
            ] = boltzmann_probability
            partition_normalization += boltzmann_probability

        for spin_states in self.boltzmann_lookup_table.keys():
            self.boltzmann_lookup_table[spin_states] /= partition_normalization

    def _find_boltzmann_probability(
            self, site_spin, number_spin_up, number_spin_down
    ):

        boltzmann_probability = self.boltzmann_lookup_table[
            (site_spin, number_spin_up, number_spin_down)
        ]

        return boltzmann_probability

    def _sweep_metropolis_slow(self, external_field, temperature):
        """Single pass of Metropolis (slow) algorithm.

        Parameters
        ----------
        external_field : float
            External magnetic field in units of (check units)

        temperature : float
            System temperature in units of (check units)

        """
        site_x, site_y, site_z = (
            random.randint(0, self.number_sites_along_xyz - 1),
            random.randint(0, self.number_sites_along_xyz - 1),
            random.randint(0, self.number_sites_along_xyz - 1),
        )
        neighbors = [
            (site_x - 1, site_y, site_z),
            (site_x + 1, site_y, site_z),
            (site_x, site_y - 1, site_z),
            (site_x, site_y - 1, site_z),
            (site_x, site_y, site_z - 1),
            (site_x, site_y, site_z + 1),
        ]

        site_spin = self[site_x, site_y, site_z]
        neighbor_spins = []
        for neighbor_x, neighbor_y, neighbor_z in neighbors:
            neighbor_spins.append(self[neighbor_x, neighbor_y, neighbor_z])

        change_in_total_energy = (
            -2.0 * site_spin * (external_field + sum(neighbor_spins))
        )

        if change_in_total_energy > temperature * math.log(random.random()):
            self[site_x, site_y, site_z] = -self[site_x, site_y, site_z]
            self.magnetization += 2 * self[site_x, site_y, site_z]

    def _write_state_snapshot_to_disk(self, external_field, temperature):
        """Save system state to output file.

        Parameters
        ----------
        external_field : float
            External magnetic field in units of (check units)

        temperature : float
            System temperature in units of (check units)

        """
        output_directory = Path(".").cwd() / "results"
        output_directory.mkdir(parents=True, exist_ok=True)
        output_filename = (
            output_directory / "final_states.csv".format(temperature)
        )

        state_full_system = []
        column_names = ["site_x", "site_y", "site_z", "spin"]
        self._gather_current_state(state_full_system)

        final_states_database = pd.DataFrame(
            data=state_full_system,
            columns=column_names,
        )
        final_states_database["temperature"] = temperature
        final_states_database["external_field"] = external_field

        _save_database_to_disk(
            database=final_states_database,
            output_filename=output_filename,
            mode="append",
            header=False if output_filename.exists() else True,
        )


def main(
        number_sites_along_xyz=10,
        steps=25000,
        external_field_sweep_start=1,
        external_field_sweep_end=11,
        temperature_sweep_start=1,
        temperature_sweep_end=11,
        algorithm="metropolis",
        seed=None,
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
        "external_field",
        "mean_magnetization",
        "std_dev_magnetization",
    ]
    output_directory = Path(".").cwd() / "results"
    output_directory.mkdir(parents=True, exist_ok=True)
    output_filename = (
        output_directory / "runs_magnetization_vs_temperature.csv"
    )
    number_simulations = (
        (external_field_sweep_end - external_field_sweep_start) *
        (temperature_sweep_end - temperature_sweep_start)
    )
    simulation_results_database = pd.DataFrame(
        columns=column_names,
        index=np.arange(number_simulations),
    )

    _external_field_and_temperature_sweep(
        ising=ising,
        steps=steps,
        column_names=column_names,
        external_field_sweep_start=external_field_sweep_start,
        external_field_sweep_end=external_field_sweep_end,
        temperature_sweep_start=temperature_sweep_start,
        temperature_sweep_end=temperature_sweep_end,
        simulation_results_database=simulation_results_database,
        algorithm=algorithm,
        seed=seed,
    )

    _save_database_to_disk(
        database=simulation_results_database,
        output_filename=output_filename,
    )


def _external_field_and_temperature_sweep(
        ising,
        steps,
        column_names,
        external_field_sweep_start,
        external_field_sweep_end,
        temperature_sweep_start,
        temperature_sweep_end,
        simulation_results_database,
        algorithm,
        seed,
):
    number_temperature_sweeps = temperature_sweep_end - temperature_sweep_start
    for loop_index, external_field in enumerate(iterable=range(
            external_field_sweep_start, external_field_sweep_end), start=0):
        simulation_results = []

        _temperature_sweep(
            ising=ising,
            steps=steps,
            external_field=external_field,
            simulation_results=simulation_results,
            sweep_start=temperature_sweep_start,
            sweep_end=temperature_sweep_end,
            algorithm=algorithm,
            seed=seed,
        )

        simulation_count_start = loop_index * number_temperature_sweeps
        simulation_count_end = (
            simulation_count_start + number_temperature_sweeps - 1
        )

        # yapf: disable
        simulation_results_database.loc[
            simulation_count_start:simulation_count_end,
            ("temperature", "mean_magnetization", "std_dev_magnetization")
        ] = simulation_results
        simulation_results_database.loc[
            simulation_count_start:simulation_count_end,
            ("external_field")
        ] = "h = {0}".format(external_field)
        # yapf: enable


def _save_database_to_disk(
        database,
        output_filename,
        mode="write",
        header=True,
):

    if mode == "write":
        database.to_csv(
            path_or_buf="{0}.gz".format(output_filename),
            index=False,
            mode="w",
            header=header,
            compression="gzip",
        )
    elif mode == "append":
        database.to_csv(
            path_or_buf="{0}.gz".format(output_filename),
            index=False,
            mode="a",
            header=header,
            compression="gzip",
        )


def _temperature_sweep(
        ising,
        steps,
        external_field,
        simulation_results,
        sweep_start,
        sweep_end,
        algorithm,
        seed,
):
    for temperature in range(sweep_start, sweep_end):
        simulation_statistics = ising.run_simulation(
            steps=steps,
            temperature=temperature,
            external_field=external_field,
            algorithm=algorithm,
            seed=seed,
        )
        simulation_results.append(simulation_statistics)


if __name__ == '__main__':

    number_sites_along_xyz = 10
    steps = 25000

    main(
        number_sites_along_xyz=number_sites_along_xyz,
        steps=steps,
    )
