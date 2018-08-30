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

g     total energy              =    -222.16750131 Ry
a     total energy              =    -222.09464546 Ry
f     total energy              =    -222.06793449 Ry
n     total energy              =    -222.06921053 Ry

MG    	 Si 	Si*Sj1 	 Si*Sj2 
NM    	0	0	0
FM    	0.5	8	6
AFM_a 	0.5	-8	6
AFM_g 	0.5	0	-2

"""

#                E_f            E_a         E_g
E = np.array([  -222.06793449 , -222.09464546  ,  -222.16750131 ]) #Ry
natoms = np.array([ 4, 4, 4])
#per atom
E = (E/natoms)# - (-54.945) 

Sj = np.array( [[ 8.0, 6.0, 1.0],       #=Efm
                [-8.0, 6.0, 1.0],       #=Eafm1
                [ 0.0,-2.0, 1.0]] )     #=Eafm2

#solve with least squares
J = np.linalg.lstsq(Sj, E)

J[0]

#convert to eV
#J = J[0] * 13.605698065894 

#J[0]
#array([  7.52373031e-03,   2.55273958e-03,  -5.52035800e+01])


