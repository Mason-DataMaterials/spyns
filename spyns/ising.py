# -*- coding: utf-8 -*-
# Copyright (c) Mason DataMaterials Group.
# Distributed under the terms of the MIT License.

"""
Created on Thu Aug 31 15:39:03 2017

@author: ssilayi
"""

import logging
import math
from random import *

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

    def __init__(self, n):

        self.n = n
        self.s = [[[1 for x in range(n)] for y in range(n)] for z in range(n)]
        self.magnetization = n**3

    def __getitem__(self, point):

        n = self.n
        x, y, z = point

        return self.s[(x + n) % n][(y + n) % n][(z + n) % n]

    def __setitem__(self, point, value):

        n = self.n
        x, y, z = point
        self.s[(x + n) % n][(y + n) % n][(z + n) % n] = value

    def configuration(self, h, t):

        n = self.n
        config = []
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    cnfg = np.array([i, j, k, self.s[i][j][k]])
                    config.append(cnfg)

        tname = "h=" + str(h) + "/configurations/t=" + str(t) + ".txt"
        np.savetxt(tname, config)

    def step(self, t, h):

        n = self.n
        x, y, z = randint(0, n - 1), randint(0, n - 1), randint(0, n - 1)
        neighbors = [(x - 1, y, z), (x + 1, y, z), (x, y - 1, z),
                     (x, y - 1, z), (x, y, z - 1), (x, y, z + 1)]
        dE = -2.0 * self[x, y, z] * \
            (h + sum(self[xn, yn, zn] for xn, yn, zn in neighbors))
        if dE > t * math.log(random()):
            self[x, y, z] = -self[x, y, z]
            self.magnetization += 2 * self[x, y, z]

        return self.magnetization


def main():

    ising = Ising(n=10)
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
            #fname = str(t)+".txt"
            #np.savetxt(fname, np.array(m))
        fname = "h=" + str(h) + "/mg.txt"
        np.savetxt(fname, data)
