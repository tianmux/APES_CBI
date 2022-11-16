#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 09:18:47 2022

@author: txin
"""
import numpy as np
import time 

from update_Vb import update_Vb_beam

def update_Vb_py(Vb, dt, Vadd, wrf, alpha ):
    for i in range(len(Vb)):
        Vb[i] = Vb[i-1]*np.exp((1j*wrf-alpha)*(dt[i]))+Vadd[i]
    return Vb

nParTot = int(1e7)
Vb1 = np.ndarray(nParTot, dtype = 'complex')
Vb2 = np.ndarray(nParTot, dtype = 'complex')
Vb1[-1] = 1
Vb2[-1] = 1

t = np.array([i for i in range(nParTot)]).astype('double')
Vadd = np.zeros(nParTot)

f = 100.0
wrf = float(2*np.pi*f)
alpha = 1e-10

start = time.time()
Vb1 = update_Vb_beam(Vb1, t, Vadd, wrf, alpha)
end = time.time()
print("Total time cython: ",end-start)
print("Time per loop cython: ",(end-start)/nParTot*1e6,"[us]")

#start = time.time()
#Vb2 = update_Vb_py(Vb2, t, Vadd, wrf, alpha)
#end = time.time()
#print("Total time pure python: ",end-start)
#print("Time per loop pure python: ",(end-start)/nParTot*1e6,"[us]")

#for i in range(len(Vb1)):
#    if Vb1[i]-Vb2[i] != 0:
#        print("Error: ", i, "    ", Vb1[i], "  ", Vb2[i])
print(Vb1[0])
print(Vb1[-1])