#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 24 01:26:53 2017
@author: Swabir
"""
import numpy as np
import matplotlib.pyplot as plt
import os,sys
#os.chdir('/Users/Swabir/Documents/Computational-Projects/autocorrelation/heisenberg/hfile/')

from scipy.optimize import curve_fit

def autocorrelation (x, x_avg) :
    """
    Compute the autocorrelation of the signal, based on the properties of the
    power spectral density of the signal.
    """
    xp = x-x_avg#np.mean(x, x_avg)
    f = np.fft.fft(xp)
    p = np.array([np.real(v)**2+np.imag(v)**2 for v in f])
    pi = np.fft.ifft(p)
    return np.real(pi)[:int(x.size/2)]/np.sum(xp**2)

def read_data(filename, m_avg):
    tmh = np.loadtxt(filename)
    t,m,h = tmh[:,0],tmh[:,1], tmh[:,2]
    m_corr = autocorrelation (m, m_avg)
    d=t[:int(len(t)/2)]
    
    return d, m_corr

  
def func(x, a, b):
    return a * np.exp(-b * x) 

def fit_tau(dir_name, data_filename):
    
    data = np.loadtxt(data_filename)
    temp, temp_avg = data[:,0], data[:,1]
    
    filename = np.sort(np.array((os.listdir(dir_name)) ).astype(int))
    tau = []
    i=0
    for fn in filename:
        if (fn in temp):
            fn = str(fn)
            
            xdata, ydata = read_data(dir_name+str(fn), temp_avg[i])
            popt, pcov = curve_fit(func, xdata, ydata)
            
            """
            plt.plot(xdata, ydata, 'b-', label='data')
            plt.plot(xdata, func(xdata, *popt), 'r-', 
                  label='fit: a=%5.3f, b=%5.3f' % tuple(popt))
            plt.show()
            
            """
            x = np.array((fn, popt[1]))
            
            tau.append(x)
            print("done " + str(fn) )
            i=i+1
        
    return tau


dir_name = "aucf_files/"
data_filename = "data.txt"
tau = fit_tau(dir_name, data_filename)
tau = np.array((tau)).astype(float)
np.savetxt("tau.txt", tau)
plt.figure()
plt.plot(tau[:,0],tau[:,1]*50,'b-', label='tau')
plt.show()
    

