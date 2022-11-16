#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 14:16:36 2022

@author: txin
"""

c_light = 299792458;
c_inver = 1/c_light;
pi = 3.141592653589793;
me = 9.1093837015e-31;
mp = 1.67262192369e-27;
qe = 1.6021766208e-19;
E0p = 938.2720813e6;
E0e = 0.510998950e6;

import pandas as pd
import numpy as np
import RF
import matplotlib.pyplot as plt
import numba as nb

qovermp = qe/mp;
qoverme = qe/me;

def Get_Inputs(fn):
    inputs = pd.read_csv(fn,sep=',')#,usecols=["Name","Values","Comments"])
    return inputs

def Gen_Pattern(inputs):
    nTrain = np.array(inputs[inputs["Name"]=="nTrain"]["Values"].iloc[0].split()).astype("int")[0]
    # the pattern is the array contains numbers is sets of two, first one is the bunches in one train, second one is the
    # number of the empty buckets, the interval is the "fill_step" in unit of bucket.
    pattern = np.array(inputs[inputs["Name"]=="Pattern"]["Values"].iloc[0].split()).astype("int")
    fill_step = np.array(inputs[inputs["Name"]=="fill_step"]["Values"].iloc[0].split()).astype("int")
    print("Number of Train: ",nTrain)
    print("Train Pattern: ",pattern)
    print("Fill step: ",fill_step)
    return nTrain, pattern, fill_step

def Init_Beam(inputs):
    nTrain,pattern,fill_step = Gen_Pattern(inputs)
    Q = np.array(inputs[inputs["Name"]=="QpB"]["Values"].iloc[0].split()).astype("float64")[0]
    
    R = np.array(inputs[inputs["Name"]=="R"]["Values"].iloc[0].split()).astype("float64")[0]
    Gamma0 = np.array(inputs[inputs["Name"]=="Gamma0"]["Values"].iloc[0].split()).astype("float64")[0]
    Ek = np.array(inputs[inputs["Name"]=="Ek"]["Values"].iloc[0].split()).astype("float64")[0]
    h  = np.array(inputs[inputs["Name"]=="h"]["Values"].iloc[0].split()).astype("float64")
    mainRF = np.array(inputs[inputs["Name"]=="mainRF"]["Values"].iloc[0].split()).astype("int")[0]
    
    Gamma = Ek/E0e
    beta = np.sqrt(1-1/Gamma**2)
    Velocity = c_light*beta
    T0 = 2*pi*R/Velocity
    f0 = 1/T0
    frf = h*f0
    Trf = 1/frf
    t0s = []
    gamma0 = []
    
    print(pattern)
    for i in range(nTrain):
        for j in range(pattern[i*2]):
            t0s.append(Trf[mainRF]/2+np.sum(pattern[:i*2])*Trf[mainRF]+j*Trf[mainRF])
    for i in range(nTrain):
        for j in range(pattern[i*2]):
            gamma0.append(Gamma)
    if 0:
        print("Charge per bunch: ",Q)
        print("Radius of the ring: ",R)
        print("Beam Ek: ",Ek/1e9, " [GeV]")
        print("Beam Gamma: ",Gamma0)
        print("Beam Ek, calculated with Ek:",Gamma)
        print("Beta of the beam: ",beta)
        print("Revolution time: ",T0," [s]")
        print("Fundamental frequency: ",f0, " [Hz]")
        print("Harmonic numbers: ",h)
        print("RF Frequencies: ",frf/1e6, " [MHz]")
        print("Time coordinates of the bunches: ", t0s)
        print("Initial Gamma of the bunches: ",gamma0)
    
    return t0s,gamma0

def TestCav(cav,N):
    fc = 1e9
    RoQ = 50
    QL = 1e3
    R = RoQ*QL

    df = 1e5
    frf = fc+df

    Ig0 = 1+0j
    Ig = Ig0
    tDriven0 = 1*1/frf
    dt = 1/frf/1e8
    cav1 = cav
    fL = cav1.wL/2/np.pi
    Z = R/(1+1j*QL*(frf/fc-fc/frf))
    Vgrec = []
    Vgref = []
    for i in range(N):
        #print("i : ",i)
        # if Ig is not changed, only thing we need to do is to propagate the phasor by rotating it.
        tDriven = (1)*tDriven0
        Ig = Ig*np.exp(1j*frf*2*np.pi*tDriven)    
        #cav1.PrintAll()

        cav1.update_Vg(Ig,tDriven,dt)
        Vgref.append(Ig*Z)
        Vgrec.append(cav1.Vg)

    tanPhi = QL*(frf/fc-fc/frf)
    Phi = np.arctan(tanPhi)
    print("Z : ",Z)
    print("TanPhi : ", QL*(frf/fc-fc/frf))
    print("Phi : ", Phi/np.pi*180)
    print("Vref : ", Vgref[-1])
    print("Vlap : ", cav1.Vg)
    print("dVg : ", cav1.Vg-Vgref[-1])
    fig,axis = plt.subplots(1,1)
    axis.plot(np.real(Vgrec[:]),'rx')
    axis.plot(np.real(Vgref[:]),'b-')
    
def ArgSort(a,b):
    ia = np.argsort(a)
    return a[ia],b[ia]

@nb.njit
def Find_idx(t,t_feed):
    iStart = 0
    iEnd = len(t)-1
    iMid = int(len(t)/2)
    nIter = 0
    if t_feed>t[iEnd]:
        return iEnd
    else:
        while iEnd-iStart>1:
            if t_feed>=t[iMid]:
                #print("in the right section.")
                iStart = iMid
            else:
                #print("in the left section.")
    
                iEnd = iMid
            iMid = iStart+int((iEnd-iStart)/2)
            nIter+=1
        #print(iStart, iEnd, iMid,nIter)
        return iStart
        
        