#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 15:39:03 2017

@author: ssilayi
"""

import numpy as np
import math
from random import *

class Ising:
    
    
    def __init__(self,n):
        self.n = n
        self.s = [[[1 for x in range(n)] for y in range(n)] for z in range (n)]
        self.magnetization = n**3
        
    def __init__(self,n):
        self.n = n
        self.s = [[[1 for x in range(n)] for y in range(n)] for z in range (n)]
        self.magnetization = n**3
        
    def __getitem__(self, point):
        n = self.n
        x,y,z = point
        return self.s[(x+n)%n][(y+n)%n][(z+n)%n]
    
    def __setitem__(self, point, value):
        n = self.n
        x,y,z = point
        self.s[(x+n)%n][(y+n)%n][(z+n)%n] = value
    
    def configuration(self, h, t):
        n = self.n
        config = [] 
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    cnfg = np.array([i, j, k , self.s[i][j][k]])
                    config.append(cnfg)
                    
        tname = "h="+str(h)+"/configurations/t="+str(t)+".txt"
        np.savetxt(tname, config)
    
    def step(self, t, h):
        n = self.n
        x,y,z = randint(0,n-1),randint(0,n-1),randint(0,n-1)
        neighbors =[(x-1,y,z),(x+1,y,z),(x,y-1,z),(x,y-1,z),(x,y,z-1),(x,y,z+1)]
        dE = -2.0*self[x,y,z]*(h+sum(self[xn,yn,zn] for xn, yn, zn in neighbors))
        if dE > t*math.log(random()):
           self[x,y,z] = -self[x,y,z]
           self.magnetization +=2*self[x,y,z]
           
        
        return self.magnetization

def main():  
    ising = Ising(n=10)
    steps = 25000
    for h in range (1,11):
        data =[]
        for t in range (1,11):
            m = [ising.step(t=float(t),h=float(h)) for k in range(steps)]
            print (h, t)
            mu = np.mean(m)
            sigma = np.std(m) 
            ising.configuration(t=t,h=h)
            x= (t, mu, sigma)
            data.append(x)    
            #fname = str(t)+".txt"    
            #np.savetxt(fname, np.array(m))
        fname = "h="+str(h)+"/mg.txt"   
        np.savetxt(fname, data)       


