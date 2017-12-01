#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 20:30:13 2017

@author: Swabir
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as so

"""
NM   1     -54.94490706 Ry
FM   1     -55.12807369 Ry
AFM1 2    -110.49690675 Ry
AFM2 4    -220.84368080 Ry


MG    	 Si 	Si*Sj1 	 Si*Sj2 
NM    	0	0	0
FM    	0.5	4	3
AFM_1 	0.5	-4	3
AFM_2 	0.5	0	-3

"""

#                Efm            Eafm1       Eafm2
E = np.array([ -55.12807369,-110.49690675,-220.87558563 ]) #Ry
natoms = np.array([ 1, 2, 4])
#per atom
E = (E/natoms)# - (-54.945) 

Sj = np.array( [[ 4.0, 3.0, 1.0],       #=Efm
                [-4.0, 3.0, 1.0],       #=Eafm1
                [ 0.0,-3.0, 1.0]] )     #=Eafm2

#solve with least squares
J = np.linalg.lstsq(Sj, E)

#convert to eV
J = J[0] * 13.605698065894 

J
#array([  2.04731206e-01,   6.94636080e-02,  -7.51083241e+02])


