#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 19:54:10 2017

@author: Swabir
"""

import numpy as np
import math
from random import *

class Heisenberg:
    
    def __init__(self,n):
        
        self.n = n
        
        self.theta = [[[1.0/np.cos(1.0-2.0*random()) for x in range(n)] for y in range(n)] for z in range (n)]
        self.phi = [[[np.pi*2.0*random() for x in range(n)] for y in range(n)] for z in range (n)]
        
        #(Sx, Sy, Sz) = (sin(theta)cos(phi), sin(theta)sin(phi), cos(theta))
        self.sx = [[[np.sin(self.theta[x][y][z])*np.cos(self.phi[x][y][z]) for x in range(n)] for y in range(n)] for z in range (n)]
        self.sy = [[[np.sin(self.theta[x][y][z])*np.sin(self.phi[x][y][z]) for x in range(n)] for y in range(n)] for z in range (n)]
        self.sz = [[[np.cos(self.theta[x][y][z]) for x in range(n)] for y in range(n)] for z in range (n)]
        self.magnetization = n**3
       
    def configuration(self, h, t):
        n = self.n
        config = [] 
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    m = np.sqrt(self.sx[i][j][k]**2. + self.sy[i][j][k]**2. + self.sz[i][j][k]**2.)
                    cnfg = np.array([i, j, k , m, self.sx[i][j][k], self.sy[i][j][k], self.sz[i][j][k]])
                    config.append(cnfg)
                    
        tname = "h="+str(h)+"/configurations/t="+str(t)+".txt"
        np.savetxt(tname, config)

    def step(self, t, h):
        
        n = self.n
        x,y,z = randint(0,n-1),randint(0,n-1),randint(0,n-1)
        neighbors =[((x+n-1)%n,y,z),((x+1)%n,y,z),(x,(y+n-1)%n,z),(x,(y+1)%n,z),(x,y,(z+n-1)%n),(x,y,(z+1)%n)]
       
        sj_x = np.sum(self.sx[xn][yn][zn] for xn, yn, zn in neighbors )
        sj_y = np.sum(self.sy[xn][yn][zn] for xn, yn, zn in neighbors )
        sj_z = np.sum(self.sz[xn][yn][zn] for xn, yn, zn in neighbors )
        
        #update theta
        dE = -2.0*np.sum([self.sx[x][y][z]*(h+sj_x), self.sy[x][y][z]*(h+sj_y), self.sz[x][y][z]*(h+sj_z)])
        
        if dE > t*math.log(random()):
            self.theta[x][y][z] = 1.0/np.cos(1.0-2.0*random())
            self.sx[x][y][z] = np.sin(self.theta[x][y][z])*np.cos(self.phi[x][y][z])
            self.sy[x][y][z] = np.sin(self.theta[x][y][z])*np.sin(self.phi[x][y][z])
            self.sz[x][y][z] = np.cos(self.theta[x][y][z])
            
            self.magnetization += 2. * np.sqrt(self.sx[x][y][z]**2. + self.sy[x][y][z]**2. + self.sz[x][y][z]**2.)
        
        #update phi
        dE = -2.0*np.sum([self.sx[x][y][z]*(h+sj_x), self.sy[x][y][z]*(h+sj_y), self.sz[x][y][z]*(h+sj_z)])
        
        if dE > t*math.log(random()):
            self.phi[x][y][z] = np.pi*2.0*random()
            self.sx[x][y][z] = np.sin(self.theta[x][y][z])*np.cos(self.phi[x][y][z])
            self.sy[x][y][z] = np.sin(self.theta[x][y][z])*np.sin(self.phi[x][y][z])
            self.sz[x][y][z] = np.cos(self.theta[x][y][z])
            
            self.magnetization += 2. * np.sqrt(self.sx[x][y][z]**2. + self.sy[x][y][z]**2. + self.sz[x][y][z]**2.)
            
        return self.magnetization
    
def main():  
    heisenberg = Heisenberg(n=4)
    steps = 125000
    for h in range (1, 11):
        data=[]
        for t in range (1,11):
            m = [heisenberg.step(t=float(t),h=float(h)) for k in range(steps)]
            print (h, t)
            mu = np.mean(m)
            sigma = np.std(m) 
            heisenberg.configuration(t=t,h=h)
            x= (t, mu, sigma)
            data.append(x)    
            #fname = str(t)+".txt"    
            #np.savetxt(fname, np.array(m))
        fname = "h="+str(h)+"/mg.txt"   
        np.savetxt(fname, data)       
     
main()        
    