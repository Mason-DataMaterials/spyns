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

import logging
import math
import random

import numpy as np

__author__ = "Swabir Silayi"
__copyright__ = "Copyright 2017, Mason DataMaterials Group"
__maintainer__ = "James Glasbrenner"
__email__ = "jglasbr2@gmu.edu"
__date__ = "August 31, 2017"

logger = logging.getLogger(__name__)


class Ising(object):
    """Monte Carlo simulation of the Ising model.

    Parameters
    ----------
    n : int
        Sets the number of sites in the square, three-dimensional Ising model.

    """

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
            ]
            for z in range(number_sites_along_xyz)
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

    def configuration(self, h, t):
        """Temporary docstring.

        Parameters
        ----------
        h : float
            Magnetic field

        t : float
            Temperature

        """
        n = self.number_sites_along_xyz
        config = []
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    cnfg = np.array(
                        [i, j, k, self.site_spin[i][j][k]]
                    )
                    config.append(cnfg)

        tname = "h={0}/configurations/t={1}.txt".format(str(h), str(t))
        np.savetxt(tname, config)

    def step(self, t, h):
        """Temporary docstring.

        Parameters
        ----------
        h : float
            Magnetic field

        t : float
            Temperature

        """
        n = self.number_sites_along_xyz
        x, y, z = (
            random.randint(0, n - 1),
            random.randint(0, n - 1),
            random.randint(0, n - 1)
        )
        neighbors = [
            (x - 1, y, z),
            (x + 1, y, z),
            (x, y - 1, z),
            (x, y - 1, z),
            (x, y, z - 1),
            (x, y, z + 1)
        ]
        dE = (
            -2.0 * self[x, y, z] * (
                h + sum(self[xn, yn, zn] for xn, yn, zn in neighbors)
            )
        )

        if dE > t * math.log(random.random()):

            self[x, y, z] = -self[x, y, z]
            self.magnetization += 2 * self[x, y, z]

        return self.magnetization


def main():
    """Temporary docstring."""
    ising = Ising(number_sites_along_xyz=10)
    steps = 25000

    for h in range(1, 11):

        data = []

        for t in range(1, 11):

            m = [ising.step(t=float(t), h=float(h)) for k in range(steps)]
            print(h, t)
            mu = np.mean(m)
            sigma = np.std(m)

            ising.configuration(t=t, h=h)
            x = (t, mu, sigma)
            data.append(x)

            # fname = str(t)+".txt"
            # np.savetxt(fname, np.array(m))

        fname = "h=" + str(h) + "/mg.txt"
        np.savetxt(fname, data)
