#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 24 10:32:08 2022

@author: txin
"""


import numpy as np

import particle
import utl
import RF
import OneTurn
import matplotlib.pyplot as plt
import time
import copy
from scipy.optimize import curve_fit


def func_exp_fit(x, a, b, c):
    return a*np.exp(b*x)+c

def plot_controids(beam, sampRate, iturn,Gamma0,bucketHieght): # iturn is the turn number we are looking at  
    fig,axis = plt.subplots(2,1)
    fig.set_figheight(16)
    fig.set_figwidth(30)
    axis[0].plot(beam.phi_cents[iturn][:],'r.',ms=10)
    axis[1].plot(beam.g_cents[iturn]-Gamma0,'r.',ms=10)
    axis[1].hlines(y=bucketHieght*Gamma0,xmin=0,xmax=len(beam.g_cents[iturn]),colors='b',linestyles='dashed')
    #axis[2].plot(np.sqrt(((beam.g_cents[iturn]-Gamma0)/Gamma0)**2+(beam.phi_cents[iturn]/np.pi)**2),'r.',ms=10)
    
    axis[0].set_xlabel("Bunch index",fontsize=30)
    axis[0].set_ylabel("phi centroids",fontsize=30)
    axis[1].set_xlabel("Bunch index",fontsize=30)
    axis[1].set_ylabel("Gamma centroids",fontsize=30)
    axis[0].set_title("phi centroinds along the bunch train at turn "+str(iturn*sampRate),fontsize=30)
    axis[1].set_title("Gamma centroinds along the bunch train at turn "+str(iturn*sampRate),fontsize=30)

    axis[0].tick_params(labelsize=30)
    axis[1].tick_params(labelsize=30)

    save_time = time.time()
    fig.savefig("./figs/Centroids"+str(save_time)+".png",bbox_inches='tight')

def get_mode_amp(beam,Gamma0):
    Amp = np.ndarray((len(beam.g_cents),len(beam.g_cents[0]))).astype(complex)
    for i in range(len(beam.g_cents)):
        Amp[i] = copy.deepcopy(np.fft.fft(beam.g_cents[i]-Gamma0))
    return Amp

def fit_growth_rate(Amp,mu,istart,iend,istep,b_gues):
    rng = np.array([istart+i*istep for i in range(iend-istart)])
    #print(([-np.abs(Amp.T[mu][istart])-1,-b_gues,-np.abs(Amp.T[mu][istart])+1], [np.abs(Amp.T[mu][istart])+1, b_gues, np.abs(Amp.T[mu][istart])+1]))
    popt, pcov = curve_fit(func_exp_fit, rng, \
                               Amp.T[mu][istart:iend],bounds=([-np.abs(Amp.T[mu][istart])-1,-b_gues,-np.abs(Amp.T[mu][istart])-1], [np.abs(Amp.T[mu][istart])+1, b_gues, np.abs(Amp.T[mu][istart])+1]),\
                               maxfev=200000000)
    return popt,pcov

def get_growth_rate(Amp,mus,istart,iend,istep,b_gues):
    Omegas_Im= []
    for i in range(len(mus)):
        rng = np.array([istart+i*istep for i in range(iend-istart)])
        popt, pcov = curve_fit(func_exp_fit, rng, \
                               Amp.T[mus[i]][istart:iend],bounds=([-np.abs(Amp.T[i][istart])-1,-b_gues,-np.abs(Amp.T[i][istart])-1], [np.abs(Amp.T[i][istart])+1, b_gues, np.abs(Amp.T[i][istart])+1]),\
                               maxfev=200000000)
        Omegas_Im.append(popt[1])
    return Omegas_Im

def get_fitted_curve(Amp,mu,istart,iend,istep,popt):
    rng = np.array([istart+i*istep for i in range(iend-istart)])
    return rng,func_exp_fit(rng,popt[0],popt[1],popt[2])

def plot_mode_amp(Amp,Amp_fit,rng,imode,istart,iend):
    fig,axis = plt.subplots(1,1)
    fig.set_figheight(16)
    fig.set_figwidth(30)
    axis.plot(rng, np.abs(Amp.T[imode])[istart:iend],'r.',ms=10)
    axis.plot(rng, Amp_fit, 'b-', ms=10)
    axis.set_title("Amplitude of mode "+str(imode),fontsize=30)

    axis.tick_params(labelsize=30)

    save_time = time.time()
    fig.savefig("./figs/Amplitude_"+str(imode)+"_"+str(save_time)+".png",bbox_inches='tight')
    
def write_output(fn, data,index):
    fh = open(fn,'w')
    for i in range(len(index)):
        fh.write(str(index[i])+',    '+str(data[i])+'\n')
    fh.close()
    
def plot_phi_diff(phases,pattern):
    dPhi = phases[:pattern[0]]-phases[-pattern[0]:]
    fig,axis = plt.subplots(1,1)
    fig.set_figheight(16)
    fig.set_figwidth(30)
    rng = range(pattern[0])
    axis.plot(rng, dPhi,'r.',ms=10)
    axis.set_title("Phase difference between colliding bunches",fontsize=30)

    axis.tick_params(labelsize=30)

    save_time = time.time()
    fig.savefig("./figs/dPhi"+"_"+str(save_time)+".png",bbox_inches='tight')
    
    return
def plot_z_diff(phases,pattern,lmd):
    dPhi = phases[:pattern[0]]-phases[-pattern[0]:]
    fig,axis = plt.subplots(1,1)
    fig.set_figheight(16)
    fig.set_figwidth(30)
    rng = range(pattern[0])
    axis.plot(rng, dPhi/360*lmd*1000.0,'r.',ms=10)
    axis.set_title("Position difference between colliding bunches [mm]",fontsize=30)

    axis.tick_params(labelsize=30)

    save_time = time.time()
    fig.savefig("./figs/dz"+"_"+str(save_time)+".png",bbox_inches='tight')
    
    return