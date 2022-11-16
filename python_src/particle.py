#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 13:44:46 2022

@author: txin
"""
import numpy as np
import copy
pi = np.pi
class Particle2D:
    
    def __init__(self,t,gamma,q):
        self.t = t
        self.gamma = gamma
        self.q = q
        
    def update(self,t,gamma,q):
        self.t = t
        self.gamma = gamma
        self.q = q

    def Print(self):
        print("Particle t: ",self.t," s.")
        print("Particle gamma: ",self.gamma,'.')
        print("Particle q: ",self.q, " C.")

class Bunch:
    def __init__(self,nPar,ts,gammas,Q,dim):
        self.nPar = 0
        self.Pars = []
        self.Vgs = []
        self.Vbs = []
        self.Q = 0
        self.nPar = nPar
        self.Q = Q
        q = self.Q/self.nPar
        if dim=="2D":
            for i in range(len(ts)):
                self.Pars.append(Particle2D(ts[i],gammas[i],q))
            self.Vgs = np.ndarray(self.nPar,dtype='complex')
            self.Vbs = np.ndarray(self.nPar,dtype='complex')
    
    def get_M1(self):
        ts = np.zeros(self.nPar)
        gs = np.zeros(self.nPar)
        for i in range(self.nPar):
            ts[i] = self.Pars[i].t
            gs[i] = self.Pars[i].gamma
        return np.average(ts),np.average(gs)
    def Print(self):
        for i in range(len(self.Pars)):
            self.Pars[i].Print()
        
class Beam:
    def __init__(self,nPb,nBunch,ts,gammas,qPb,dim,h0,pattern, fill_step, nSamp):
        self.nPar = 0
        self.ts = []
        self.gammas = []
        self.qs = []
        self.Vgs = []
        self.Vbs = []
        self.Vadds = []
        self.nPb = nPb
        self.nBunch = nBunch
        self.nTot = nPb*nBunch

        self.Q = qPb
        self.buckets = np.zeros(int(h0[0])) # store the status of having bunch or not. 
        self.t_cents = np.zeros((nSamp,nBunch))
        self.phi_cents = np.zeros((nSamp,nBunch))

        self.g_cents = np.zeros((nSamp,nBunch))
        
        current_bucket = 0
        for i in range(int(len(pattern)/2)):
            for j in range(pattern[i*2]):
                self.buckets[current_bucket] = 1
                current_bucket += fill_step
            current_bucket +=pattern[i*2+1]*fill_step
            
        if dim=="2D":
            self.ts = copy.deepcopy(ts)
            self.gammas = copy.deepcopy(gammas)
            self.qs = np.zeros(nPb*nBunch)
            
        for i in range(self.nTot):
            self.qs[i] = self.Q/self.nPb
            
    def get_M1(self,BN,nPar):
        return np.average(self.ts[BN*nPar:(BN+1)*nPar]),np.average(self.gammas[BN*nPar:(BN+1)*nPar])
    
    def get_M1s(self,iSmp,Trf,t,t0,tCentroids):
        for i in range(self.nBunch):
            self.t_cents[iSmp][i] = np.average(self.ts[i*self.nPb:(i+1)*self.nPb])
            self.g_cents[iSmp][i] = np.average(self.gammas[i*self.nPb:(i+1)*self.nPb])
        self.phi_cents[iSmp] = (self.t_cents[iSmp]-t-tCentroids+t0)/Trf[0]*360
        return 
    def Print(self):
        return