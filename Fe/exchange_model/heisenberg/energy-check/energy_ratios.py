#!/bin/env python3
import os, sys
import numpy as np

j13 = np.loadtxt("energies.13")
j1 = np.loadtxt("energies.1")
E_13 = np.average(j13[:,2]) 
E_1 = np.average(j1[:,2]) 
E_factor = E_13/E_1
print("Energy Ratio: " + str(E_factor) + " in Ry")
