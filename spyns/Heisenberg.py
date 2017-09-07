#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 21:03:41 2017

@author: Swabir
"""
import os, sys
import numpy as np
import matplotlib.pyplot as plt
from random import *
import math


class Heisenberg:
    
    def __init__(self, n):
        
        self.n = n
        
        self.theta = [[[1.0/np.cos(1.0-2.0*random()) for x in range(n)] for y in range(n)] for z in range (n)]
        self.phi = [[[np.pi*2.0*random() for x in range(n)] for y in range(n)] for z in range (n)]
        
        #(Sx, Sy, Sz) = (sin(theta)cos(phi), sin(theta)sin(phi), cos(theta))
        #theta, phi = self.theta, self.phi
    
        self.sx = [[[np.sin(self.theta[x][y][z])*np.cos(self.phi[x][y][z]) for x in range(n)] for y in range(n)] for z in range (n)]
        self.sy = [[[np.sin(self.theta[x][y][z])*np.sin(self.phi[x][y][z]) for x in range(n)] for y in range(n)] for z in range (n)]
        self.sz = [[[np.cos(self.theta[x][y][z]) for x in range(n)] for y in range(n)] for z in range (n)]
        
        self.M = 0.0
        self.H = 0.0
    
    def magnetization (self):
        
        return np.sum(np.square((1./self.n)*np.array((np.sum(self.sx), np.sum(self.sy), np.sum(self.sz)))))
    
    def energy (self, h):
        
        n = self.n
        H_step = 0.0
        for x in range (n):
            for y in range(n):
                for z in range(n):
                    
                    sx = np.sin(self.theta[x][y][z])*np.cos(self.phi[x][y][z]) 
                    sy = np.sin(self.theta[x][y][z])*np.sin(self.phi[x][y][z]) 
                    sz = np.cos(self.theta[x][y][z]) 
                    
                    Si = np.array((sx, sy, sz))
                    
                    neighbors =[((x-1)%n,y,z),((x+1)%n,y,z),(x,(y-1)%n,z),(x,(y-1)%n,z),(x,y,(z-1)%n),(x,y,(z+1)%n)]
                    sxn = sum(self.sx[xn][yn][zn] for xn, yn, zn in neighbors)
                    syn = sum(self.sy[xn][yn][zn] for xn, yn, zn in neighbors)
                    szn = sum(self.sz[xn][yn][zn] for xn, yn, zn in neighbors)
                    
                    Sj = np.array((sxn, syn, szn))
                    
                    H_step += -np.dot(Si , h+Sj)
        
        return 0.5*H_step
        
        
        
        
    def step(self, t, h):
        
        n = self.n
        
        #pick a point
        x,y,z = randint(0,n-1),randint(0,n-1),randint(0,n-1)
        
        #find neighbors
        neighbors =[((x-1)%n,y,z),((x+1)%n,y,z),(x,(y-1)%n,z),(x,(y-1)%n,z),(x,y,(z-1)%n),(x,y,(z+1)%n)]
        
        
        #update theta holding phi constant
        #possible new theta
        theta = 1.0/np.cos(1.0-2.0*random())
        
        sx = np.sin(theta)*np.cos(self.phi[x][y][z]) 
        sy = np.sin(theta)*np.sin(self.phi[x][y][z]) 
        sz = np.cos(theta) 
        
        sxn = sum(self.sx[xn][yn][zn] for xn, yn, zn in neighbors)
        syn = sum(self.sy[xn][yn][zn] for xn, yn, zn in neighbors)
        szn = sum(self.sz[xn][yn][zn] for xn, yn, zn in neighbors)
        
        #deltaE = -J (Si_new – Si_old)*(S1 + S2+ S3+ S4+ S5+ S6)
        Si_new = np.array((sx, sy, sz))
        Si_old = np.array((self.sx[x][y][z], self.sy[x][y][z],self.sz[x][y][z] ))
        Sj = np.array((sxn, syn, szn))
        
        deltaE = -2.* np.dot((Si_new - Si_old), h+Sj)
        if deltaE > t*math.log(random()):
            self.theta [x][y][z]= theta
            
        #update phi holding theta constant
        #possible new phi 
        phi = np.pi*2.0*random() 
        sx = np.sin(self.theta[x][y][z])*np.cos(phi) 
        sy = np.sin(self.theta[x][y][z])*np.sin(phi) 
        sz = np.cos(self.theta[x][y][z]) 
        
        sxn = sum(self.sx[xn][yn][zn] for xn, yn, zn in neighbors)
        syn = sum(self.sy[xn][yn][zn] for xn, yn, zn in neighbors)
        szn = sum(self.sz[xn][yn][zn] for xn, yn, zn in neighbors)
        
        #deltaE = -J (Si_new – Si_old)*(S1 + S2+ S3+ S4+ S5+ S6)
        Si_new = np.array((sx, sy, sz))
        Si_old = np.array((self.sx[x][y][z], self.sy[x][y][z],self.sz[x][y][z] ))
        
        Sj = np.array((sxn, syn, szn))
        
        deltaE = -2.* np.dot((Si_new - Si_old), h+Sj)
        
        if deltaE > t*math.log(random()):
            self.phi [x][y][z]= phi
         
            
        #update configuration
        self.sx = [[[np.sin(self.theta[x][y][z])*np.cos(self.phi[x][y][z]) for x in range(n)] for y in range(n)] for z in range (n)]
        self.sy = [[[np.sin(self.theta[x][y][z])*np.sin(self.phi[x][y][z]) for x in range(n)] for y in range(n)] for z in range (n)]
        self.sz = [[[np.cos(self.theta[x][y][z]) for x in range(n)] for y in range(n)] for z in range (n)]
        
        #energy and magnetization after each step
        self.M = self.magnetization()
        self.H = self.energy(h)
        
        MH = self.M, self.H
        
        return MH
            
        
        
Heiss = Heisenberg(4)

mcSteps = 3000

h = 1.0
T = np.linspace(1,10, 12)
x = []
for t in (T):
    MH = [Heiss.step(t,h) for k in range(mcSteps)]
    mh = np.array((MH))
    m = np.mean(mh[:, 0])
    h = np.mean(mh[:, 1])
    x.append((t, m, h))
    
#plt.figure
#plt.plot(np.array(x[:, 0]), label="Mi")
#plt.plot(np.array(x[:, 1]), label="Hi")    
fname = str(h)+".txt"
np.savetxt(fname, np.array(x))

    
    
    
 
    