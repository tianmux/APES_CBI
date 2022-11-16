#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 13:09:06 2022

@author: txin
"""

import numpy as np
import matplotlib.pyplot as plt

def getZ(w0,w,R,Q):
    return R*w/(w+1j*Q*(w0-w**2/w0))

pi = np.pi
f0 = 1e9
n_w_oneside = 50000
df = 1e0

h = 100000
frev = f0/h
detune = 1e4
w0 = (f0-detune)*2*pi

Trev = 1/frev
f = np.array([(f0+frev)*1+df*(i-n_w_oneside) for i in range(n_w_oneside*2+1)])
fm = np.array([-f0-frev+df*(i-n_w_oneside) for i in range(n_w_oneside*2+1)])

w = f*2*pi
wm = fm*2*pi
dt = -Trev*0.0
T = Trev+dt
R = 1e9
Q = 4e4

Z = getZ(w0,w,R,Q)
Zm = getZ(w0,wm,R,Q)
epsilon = 0
gc = 1/R*1
alpha = gc*(1-epsilon/(1-epsilon))

print(alpha)
Z_comb = Z/(1-Z*alpha*np.exp(-1j*w*T))
Zm_comb = Zm/(1-Zm*alpha*np.exp(-1j*wm*T))

fig,axis = plt.subplots(1,1)
axis.plot(f,np.log10(np.real(Z)),'r.-')
axis.plot(f,np.log10(np.real(Z_comb)),'b.-')

#axis.plot(f,np.log10(np.real(Zm)),'g.-')
#axis.plot(f,np.log10(np.real(Zm_comb)),'y.-')

#axis.axvline(x = f0+frev*0.5)
#axis.axvline(x = f0+frev*1.5)

#axis.axvline(x = 0)
#axis.axvline(x = 0+1/T)
#axis.axvline(x = 0+1/T*2)

axis.axvline(x = f0,color='r')
axis.axvline(x = f0+frev,color='r')

axis.axvline(x = f0+1/T,color='b')
axis.axvline(x = f0+1/T*2,color='b')
axis.set_ylim(5)
fig.set_figheight(10*1)
fig.set_figwidth(10)
axis.tick_params(labelsize=10)
