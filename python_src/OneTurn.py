#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 16:07:50 2022

@author: txin
"""
import numpy as np
import utl

def oneTurnMap_beam(GMTSQ,Gamma0,gamma,t,T0,dynamic_on):
    eta = 1/GMTSQ-1/Gamma0**2
    delta = gamma/Gamma0-1
    t += T0*(1+delta*eta*dynamic_on)

def rad_map_beam(vRad,Gamma0, gamma,rad_on,damp_coeff,excit_coeff):
    tmp_rand = np.random.uniform(-3,1,len(gamma))
    if rad_on:         
        gamma -= vRad[0]/utl.E0e
        gamma -= (gamma-Gamma0)*damp_coeff
        gamma -= excit_coeff*(tmp_rand)
def oneTurnMap(GMTSQ,Gamma0,bunch,T0,dynamic_on):
    
    eta = 1/GMTSQ-1/Gamma0**2
    gammas = np.ndarray(len(bunch.Pars),dtype="float")
    ts = np.ndarray(len(bunch.Pars),dtype="float")
    for i in range(len(bunch.Pars)):
        gammas[i] = bunch.Pars[i].gamma
        ts[i] = bunch.Pars[i].t
    delta = gammas/Gamma0-1
    ts += T0*(1+delta*eta)
    for i in range(len(bunch.Pars)):
        bunch.Pars[i].t = ts[i]
    
def rad_map(vRad,Gamma0, bunch,rad_on,damp_coeff):
    if rad_on:         
        for i in range(bunch.nPar):
            bunch.Pars[i].gamma-=vRad/utl.E0e
            bunch.Pars[i].gamma-=(bunch.Pars[i].gamma-Gamma0)*damp_coeff