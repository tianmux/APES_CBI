#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 11:59:36 2022

@author: txin
"""

import numpy as np

import particle
import utl
import RF
import OneTurn
import matplotlib.pyplot as plt
import time


#==============================================================================

t0 = 0

inputs = utl.Get_Inputs('input_Higgs_Booster.txt')
utl.Init_Beam(inputs)
t0s = utl.Init_Beam(inputs)
q = np.array(inputs[inputs["Name"]=="QpB"]["Values"].iloc[0].split()).astype("float")[0]
sigt = np.array(inputs[inputs["Name"]=="sigt"]["Values"].iloc[0].split()).astype("float")[0]
Gamma0 = np.array(inputs[inputs["Name"]=="Gamma0"]["Values"].iloc[0].split()).astype("float")[0]
siggamma = np.array(inputs[inputs["Name"]=="siggamma"]["Values"].iloc[0].split()).astype("float")[0]
nPar = np.array(inputs[inputs["Name"]=="Npar"]["Values"].iloc[0].split()).astype("int")[0]
qPb = np.array(inputs[inputs["Name"]=="QpB"]["Values"].iloc[0].split()).astype("float")[0]
nBunch = np.array(inputs[inputs["Name"]=="nBunch"]["Values"].iloc[0].split()).astype("int")[0]
fillStep = np.array(inputs[inputs["Name"]=="fill_step"]["Values"].iloc[0].split()).astype("int")[0]

rRing = np.array(inputs[inputs["Name"]=="R"]["Values"].iloc[0].split()).astype("float")[0]
h = np.array(inputs[inputs["Name"]=="h"]["Values"].iloc[0].split()).astype("int")[0]
nStation =  np.array(inputs[inputs["Name"]=="nStation"]["Values"].iloc[0].split()).astype("int")[0]
RoQ = np.array(inputs[inputs["Name"]=="RoQ"]["Values"].iloc[0].split()).astype("float")[0]
RoQ = RoQ/nStation

QL = np.array(inputs[inputs["Name"]=="QL"]["Values"].iloc[0].split()).astype("float")[0] 
R = RoQ*QL
vRad = np.array(inputs[inputs["Name"]=="vRad"]["Values"].iloc[0].split()).astype("float")[0]
vQuad =  np.array(inputs[inputs["Name"]=="vQuad"]["Values"].iloc[0].split()).astype("float")[0] 
vRad = vRad/nStation
vSync = vRad
vQuad = vQuad/nStation

df =  np.array(inputs[inputs["Name"]=="df"]["Values"].iloc[0].split()).astype("float")[0] 
nWatch = np.array(inputs[inputs["Name"]=="nWatch"]["Values"].iloc[0].split()).astype("int")[0] 
nSamp = np.array(inputs[inputs["Name"]=="nSamp"]["Values"].iloc[0].split()).astype("int")[0] 

T0 = 2*np.pi*rRing/utl.c_light
frf = 1/T0*h

fc = frf+df
Trf = 1/frf
tDriven0 = 1*1/frf*h

cav1 = RF.Cav(fc,frf,RoQ,QL,0)
cav2 = RF.Cav(fc,frf,RoQ,QL,0)

fL = cav1.wL/2/np.pi
Z = R/(1+1j*QL*(frf/fc-fc/frf))
Vc = vSync+1j*vQuad
IbDC = -qPb*nBunch/T0
Ig = Vc/Z-IbDC*2

cav1.Igm1 = Ig
cav1.Vgm1 = Ig*Z
cav1.Ugm1 = Ig*Z/1j/cav1.wrf
cav1.update_Ig(Ig,t0)

cav2.Igm1 = Ig
cav2.Vgm1 = Ig*Z
cav2.Ugm1 = Ig*Z/1j/cav2.wrf
cav2.update_Ig(Ig,t0)
n_wait = int(1e9)
t0 = (n_wait)*Trf

Ib = qPb*nBunch/T0

tCentroids = np.array([t0+Trf*i*fillStep for i in range(nBunch)])


gammas = np.array(np.random.normal(Gamma0,siggamma,nPar))
dim = "2D"
#bunch1 = particle.Bunch(nPar,ts, gammas, qPb,dim)
bunches = []

for i in range(nBunch):
    #ts = np.array(np.random.normal(tCentroids[i],sigt,nPar))
    ts = np.array([tCentroids[i]+0*Trf/1e2 for _ in range(nPar)])
    ts = np.sort(ts)
    bunches.append(particle.Bunch(nPar,ts,gammas,qPb,dim))
    #print(len(bunches[-1].Pars))
#==============================================================================
t_start = time.time()
debug0 = 0
if debug0:
    N = int(1e5)
    Vgrec = []
    Vgref = []
    for i in range(N):
        
        dt = i*tDriven0
        #cav1.PrintAll()
        cav1.update_Vg(dt)
        #cav2.update_Vg(dt)

        Vgref.append(Ig*Z*(1-np.exp(-cav1.alpha*dt)))

        Vgrec.append(cav1.calcul_Vg(dt))
    
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

#==============================================================================

# Multi-turn kick
tWatch = np.ndarray(shape=(nWatch,nSamp,nBunch),dtype='float')
gammaWatch = np.ndarray(shape=(nWatch,nSamp,nBunch),dtype='float')

debug2 = 1
if debug2:
    nRamp = int(400)
    tRec = np.zeros(nRamp)
    gRec = np.zeros(nRamp)
    tBunch = np.zeros(nBunch)
    
    dynamic_on = 0     
    rad_on = 0
    t = t0
    for i in range(nRamp):
        print("i = ",i)
        
        cav1.update_Vg(t)
        cav2.update_Vg(t)
        for j in range(nBunch):
            cav1.kick_par(bunches[j],dynamic_on)
            OneTurn.rad_map(vRad/nStation,Gamma0, bunches[j],rad_on)
            OneTurn.oneTurnMap(inputs, bunches[j], T0/nStation)
            if nStation == 2:
                cav2.kick_par(bunches[j],dynamic_on)
                OneTurn.rad_map(vRad/nStation,Gamma0, bunches[j],rad_on)
                OneTurn.oneTurnMap(inputs, bunches[j], T0/nStation)
        t = t+T0
        tRec[i],gRec[i] = bunches[0].get_M1()
        tRec[i] -= t-t0+tCentroids[0]
        gRec[i] -=Gamma0
        cav1.update_Ig(Ig,t)
    fig,axis = plt.subplots(1,1)
    axis.plot(tRec[:],gRec[:],'r.',ms=1)
    
    dynamic_on = 1
    rad_on = 1
    nTurn = int(2000)
    tRec1 = np.zeros(nTurn)
    gRec1 = np.zeros(nTurn)
    Vq = -cav1.RoQ*cav1.wL*qPb/bunches[j].nPar
    if dynamic_on and rad_on:
        for i in range(nTurn):
             print("i = ",i)
             cav1.update_Vg(t)
             for j in range(nBunch):
                 cav1.kick_par(bunches[j],dynamic_on)
                 OneTurn.rad_map(vRad,Gamma0,bunches[j],rad_on)
                 OneTurn.oneTurnMap(inputs, bunches[j], T0/nStation)
                 if nStation == 2:
                     cav2.kick_par(bunches[j],dynamic_on)
                     OneTurn.rad_map(vRad,Gamma0, bunches[j],rad_on)
                     OneTurn.oneTurnMap(inputs, bunches[j], T0/nStation)
                 
             t = t+T0
             tRec1[i],gRec1[i] = bunches[0].get_M1()
             tRec1[i] -= t-t0+tCentroids[0]
             cav1.update_Ig(Ig,t)
             
    gRec1-=Gamma0
    fig,axis = plt.subplots(1,1)
    axis.plot(tRec1[:],gRec1[:],'r.-',ms=1)
    axis.plot(tRec1[:10],gRec1[:10],'g.-',ms=5)
    axis.plot(tRec1[-10:],gRec1[-10:],'b.-',ms=5)
    for i in range(nBunch):
        tBunch[i] = bunches[i].get_M1()[0]-t-tCentroids[i]+t0
    phiBunch = tBunch/Trf*360
    
    fig,axis = plt.subplots(1,1)
    fig.set_figheight(16)
    fig.set_figwidth(30)
    axis.plot(phiBunch,'r.',ms=10)
    axis.set_xlabel("Bunch number",fontsize=30)
    axis.set_ylabel("Phase (degree)",fontsize=30)
    axis.set_title("Phase shift along the train in first RF station.",fontsize=30)
    axis.tick_params(labelsize=30)
    save_time = time.time()
    fig.savefig("Phase_Shift_"+str(save_time)+".png",bbox_inches='tight')
    
    print("Total time : ", time.time()-t_start)
    print("Time/step : ", (time.time()-t_start)/nTurn)
    