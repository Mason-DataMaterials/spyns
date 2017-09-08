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

__author__ = "Swabir Silayi"
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

    def configuration(self, external_field, temperature):
        """Save system state to output file.

        Parameters
        ----------
        external_field : float
            External magnetic field in units of (check units)

        temperature : float
            System temperature in units of (check units)

        """
        number_sites_along_xyz = self.number_sites_along_xyz
        state_full_system = []
        for i in range(number_sites_along_xyz):
            for j in range(number_sites_along_xyz):
                for k in range(number_sites_along_xyz):
                    state_single_site = np.array(
                        [i, j, k, self.site_spin[i][j][k]]
                    )
                    state_full_system.append(state_single_site)

        state_full_system = pd.DataFrame(
            data=np.array(state_full_system),
            columns=["site_x", "site_y", "site_z", "spin"]
        )

        output_directory = (
            Path(".").cwd() /
            "results/h_{0}/configurations".format(external_field)
        )
        output_directory.mkdir(parents=True, exist_ok=True)
        output_filename = output_directory / "t_{0}.csv".format(temperature)

        state_full_system.to_csv(
            path_or_buf="{0}".format(output_filename),
            index=False,
        )

    def step(self, external_field, temperature):
        """Advance the Monte Carlo simulation by one time step.

        Parameters
        ----------
        external_field : float
            External magnetic field in units of (check units)

        temperature : float
            System temperature in units of (check units)

        """
        n = self.number_sites_along_xyz
        site_x, site_y, site_z = (
            random.randint(0, n - 1),
            random.randint(0, n - 1),
            random.randint(0, n - 1),
        )
        neighbors = [
            (site_x - 1, site_y, site_z),
            (site_x + 1, site_y, site_z),
            (site_x, site_y - 1, site_z),
            (site_x, site_y - 1, site_z),
            (site_x, site_y, site_z - 1),
            (site_x, site_y, site_z + 1),
        ]
        change_in_total_energy = (
            -2.0 * self[site_x, site_y, site_z] * (
                external_field + sum(
                    self[neighbor_x, neighbor_y, neighbor_z]
                    for neighbor_x, neighbor_y, neighbor_z in neighbors
                )
            )
        )

        if change_in_total_energy > temperature * math.log(random.random()):

            self[site_x, site_y, site_z] = -self[site_x, site_y, site_z]
            self.magnetization += 2 * self[site_x, site_y, site_z]

        return self.magnetization


def main(
        number_sites_along_xyz=10,
        steps=25000,
        external_field_sweep_start=1,
        external_field_sweep_end=11,
        temperature_sweep_start=1,
        temperature_sweep_end=11
):
    """Run Ising model simulation."""
    ising = Ising(number_sites_along_xyz=number_sites_along_xyz)

    for external_field in range(external_field_sweep_start,
                                external_field_sweep_end):

        data = []

        for temperature in range(temperature_sweep_start,
                                 temperature_sweep_end):

            magnetization_history = [
                ising.step(
                    temperature=float(temperature),
                    external_field=float(external_field)
                ) for k in range(steps)
            ]
            logger.info(external_field, temperature)
            mean_magnetization = np.mean(magnetization_history)
            std_dev_magnetization = np.std(magnetization_history)

            ising.configuration(
                temperature=temperature, external_field=external_field
            )
            simulation_statistics = (
                temperature, mean_magnetization, std_dev_magnetization
            )
            data.append(simulation_statistics)

            # fname = str(t)+".txt"
            # np.savetxt(fname, np.array(m))

        data = pd.DataFrame(
            data=np.array(data),
            columns=[
                "temperature",
                "mean_magnetization",
                "std_dev_magnetization",
            ]
        )

        output_directory = (
            Path(".").cwd() / "results/h_{0}".format(external_field)
        )
        output_directory.mkdir(parents=True, exist_ok=True)
        output_filename = (output_directory / "mg.csv")

        data.to_csv(path_or_buf="{0}".format(output_filename), index=False)


if __name__ == '__main__':

    number_sites_along_xyz = 10
    steps = 25000

    main(
        number_sites_along_xyz=number_sites_along_xyz,
        steps=steps,
    )
