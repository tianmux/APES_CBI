#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 16:59:26 2022

@author: tianmu
"""


import numpy as np

import matplotlib.pyplot as plt
import copy
from scipy.optimize import curve_fit
def func_exp_fit(x, a, b, c):
    return a*np.exp(b*x)+c

def fit_growth_rate(Amp,istart,iend,istep,b_gues):
    rng = np.array([istart+i*istep for i in range(iend-istart)])
    popt, pcov = curve_fit(func_exp_fit, rng, \
                               Amp[istart:iend],bounds=([-np.abs(Amp[istart])-1,-b_gues,-np.abs(Amp[istart])-1], [np.abs(Amp[istart])+1, b_gues, np.abs(Amp[istart])+1]),\
                               maxfev=200000000)
    return popt,pcov

pi = np.pi

T0 = 1e-3
f0 = 1/T0
N_par = 10000 # number of bunches
mu = 1
fs = f0*mu
dt = T0/N_par
N_turn = 20
N_samp = N_par*N_turn

n_mu = 6
mu_offset = 1
a = np.zeros(n_mu)
for i in range(n_mu):
    a[i] = np.random.normal(0.03,0.02)
ws = fs*2*pi*np.array([i+mu_offset for i in range(n_mu)])
rand_error = 0.0

t = np.array([i*dt for i in range(N_samp)])

phi0 = 1.2*pi

def func1(a, ws, n_mu,t,turn_idx,phi0):
    tmp = np.zeros(len(t))
    for i in range(n_mu):
        tmp+=np.sin(ws[i]*t+phi0)*np.exp(a[i]*turn_idx)
    return tmp+np.random.uniform(-rand_error,rand_error,len(t))
def func2(a, ws, n_mu,t,turn_idx,phi0):
    tmp = np.zeros(len(t))
    for i in range(n_mu):
        tmp+=np.cos(ws[i]*t+phi0)*np.exp(a[i]*turn_idx)
    return tmp+np.random.uniform(-rand_error,rand_error,len(t))

turn_idx = np.array(t/T0).astype(int)
y1 = func1(a,ws,n_mu,t,turn_idx,phi0)
y2 = func2(a,ws,n_mu,t,turn_idx,phi0)

#FFT of the data 
Fy1 = np.ndarray((N_turn,N_par)).astype(complex)
for i in range(N_turn):
    Fy1[i] = copy.deepcopy(np.fft.fft(y1[i*N_par:(i+1)*N_par]))
    
# Mode info
My1 = Fy1.T

# fit the mode info
istart = 1
iend = 20
istep = 1
b_gues = 100
n_mu_fit = 10
a_fit1 = np.zeros(n_mu_fit)
a_fit2 = np.zeros(n_mu_fit)

rng = np.array(range(istart,iend,1))
y_fit1 = np.zeros((n_mu_fit,iend-istart))
y_fit2 = np.zeros((n_mu_fit,iend-istart))
for i in range(n_mu_fit):
    popt,pcov = fit_growth_rate(np.imag(My1[i+mu_offset]), istart, iend, istep, b_gues)
    a_fit1[i] = popt[1]
    y_fit1[i] = func_exp_fit(rng, popt[0], popt[1], popt[2])
    
    popt,pcov = fit_growth_rate(np.real(My1[i+mu_offset]), istart, iend, istep, b_gues)
    a_fit2[i] = popt[1]
    y_fit2[i] = func_exp_fit(rng, popt[0], popt[1], popt[2])
# Plot the original signal
rng1 = 10
rng2 = 12
fig,axis = plt.subplots(2,1)
fig.set_figheight(20)
fig.set_figwidth(10)
axis[0].plot(t[rng1*N_par:rng2*N_par], y1[rng1*N_par:rng2*N_par],'r.',ms=10)
axis[0].plot(t[rng1*N_par:rng2*N_par], y2[rng1*N_par:rng2*N_par],'b.',ms=10)
axis[1].plot(y1[rng1*N_par:rng2*N_par],y2[rng1*N_par:rng2*N_par],'g.',ms = 10)
axis[0].set_title("Original Signal.",fontsize=30)
axis[0].tick_params(labelsize=30)
axis[1].tick_params(labelsize=30)

# Plot the FFT signal
rng1 = 0
rng2 = 30
fig,axis = plt.subplots(1,1)
fig.set_figheight(16)
fig.set_figwidth(30)
axis.plot(np.abs(np.real(Fy1[N_turn-1])[rng1:rng2]),'r.',ms=10)
axis.plot(np.abs(np.imag(Fy1[N_turn-1])[rng1:rng2]),'b.',ms=10)
axis.set_title("FFT Signal.",fontsize=30)
axis.tick_params(labelsize=30)

# Plot the Mode info
rng1 = 0
rng2 = N_turn
mu_plot = 30
fig,axis = plt.subplots(mu_plot,1)
fig.set_figheight(10*mu_plot)
fig.set_figwidth(10)
for i in range(mu_plot):
    axis[i].plot(np.abs(np.real(My1[i])[rng1:rng2]),'r.',ms=10)
    axis[i].plot(np.abs(np.imag(My1[i])[rng1:rng2]),'b.',ms=10)
    #axis[i].tick_params(labelsize=30)
axis[0].set_title("FFT Signal.",fontsize=30)


# Plot the Mode info + fitted curve
rng1 = 0
rng2 = N_turn
fig,axis = plt.subplots(n_mu_fit,1)
fig.set_figheight(10*n_mu_fit)
fig.set_figwidth(10)
for i in range(n_mu_fit):
    axis[i].plot(rng,np.abs(np.real(My1[i+mu_offset])[istart:iend]),'r.',ms=10)
    axis[i].plot(rng,np.abs(np.imag(My1[i+mu_offset])[istart:iend]),'b.',ms=10)
    axis[i].plot(rng,np.abs(y_fit1[i]),'g-',ms=10)
    axis[i].plot(rng,np.abs(y_fit2[i]),'y-',ms=10)

    #axis[i].tick_params(labelsize=30)
axis[0].set_title("FFT Signal.",fontsize=30)
print("{0:^20}{1:^20}{2:^20}{3:^20}{4:^20}".format("a","fitted a1","fitted a2","relative error1","relative error2"))
for i in range(len(a)):
    print("{0:^20.4e}{1:^20.4e}{2:^20.4e}{3:^20.4e}{4:^20.4e}".format(a[i],a_fit1[i],a_fit2[i],(a[i]-a_fit1[:len(a)][i])/a[i],(a[i]-a_fit2[:len(a)][i])/a[i]))
