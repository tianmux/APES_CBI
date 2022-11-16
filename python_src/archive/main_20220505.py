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

inputs = utl.Get_Inputs('input.txt')
utl.Init_Beam(inputs)
t0s = utl.Init_Beam(inputs)
q = np.array(inputs[inputs["Name"]=="QpB"]["Values"].iloc[0].split()).astype("float")[0]
sigt = np.array(inputs[inputs["Name"]=="sigt"]["Values"].iloc[0].split()).astype("float")[0]
Gamma0 = np.array(inputs[inputs["Name"]=="Gamma0"]["Values"].iloc[0].split()).astype("float")[0]
GMTSQ = np.array(inputs[inputs["Name"]=="GMTSQ"]["Values"].iloc[0].split()).astype("float")[0]
siggamma = np.array(inputs[inputs["Name"]=="siggamma"]["Values"].iloc[0].split()).astype("float")[0]
nPar = np.array(inputs[inputs["Name"]=="Npar"]["Values"].iloc[0].split()).astype("int")[0]
qPb = np.array(inputs[inputs["Name"]=="QpB"]["Values"].iloc[0].split()).astype("float")[0]
nTrain = np.array(inputs[inputs["Name"]=="nTrain"]["Values"].iloc[0].split()).astype("int")[0]
pattern = np.array(inputs[inputs["Name"]=="Pattern"]["Values"].iloc[0].split()).astype("int")

nBunch = np.array(inputs[inputs["Name"]=="nBunch"]["Values"].iloc[0].split()).astype("int")[0]
fillStep = np.array(inputs[inputs["Name"]=="fill_step"]["Values"].iloc[0].split()).astype("int")[0]

rRing = np.array(inputs[inputs["Name"]=="R"]["Values"].iloc[0].split()).astype("float")[0]


nStation =  np.array(inputs[inputs["Name"]=="nStation"]["Values"].iloc[0].split()).astype("int")[0]
h = np.array(inputs[inputs["Name"]=="h"]["Values"].iloc[0].split()).astype("int")[0]
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

cavs = []
for i in range(nStation):
    cavs.append(RF.Cav(fc,frf,RoQ,QL,0))

cav1 = RF.Cav(fc,frf,RoQ,QL,0)
cav2 = RF.Cav(fc,frf,RoQ,QL,0)

fL = cavs[0].wL/2/np.pi
Z = R/(1+1j*QL*(frf/fc-fc/frf))
Vc = vSync+1j*vQuad
IbDC = -qPb*nBunch/T0
Ig = Vc/Z-IbDC*2

for i in range(nStation):
    cavs[i].Igm1 = Ig
    cavs[i].Vgm1 = Ig*Z
    cavs[i].Ugm1 = Ig*Z/1j/cavs[i].wrf
    cavs[i].update_Ig(Ig,t0)
    
cav1.Igm1 = Ig
cav1.Vgm1 = Ig*Z
cav1.Ugm1 = Ig*Z/1j/cav1.wrf
cav1.update_Ig(Ig,t0)

cav2.Igm1 = Ig
cav2.Vgm1 = Ig*Z
cav2.Ugm1 = Ig*Z/1j/cav2.wrf
cav2.update_Ig(Ig,t0)

n_wait = int(1e9)
t0 = n_wait*Trf

Ib = qPb*nBunch/T0

tCentroids = np.zeros(nBunch)
ibunch = 0
ibucket = 0
for i in range(nTrain):
    for j in range(pattern[i*2]):
        tCentroids[ibunch] = t0+Trf*ibucket 
        ibunch += 1
        ibucket += fillStep
    ibucket += fillStep*pattern[i*2+1]

gammas = np.array(np.random.normal(Gamma0,siggamma,nPar))
dim = "2D"
#bunch1 = particle.Bunch(nPar,ts, gammas, qPb,dim)
bunches = []
ts_beam = np.zeros(nBunch*nPar)
gs_beam = np.zeros(nBunch*nPar)
for i in range(nBunch):
    ts = np.array(np.random.normal(tCentroids[i],sigt,nPar))
    #ts = np.array([tCentroids[i]+0*Trf/1e2 for _ in range(nPar)])
    ts = np.sort(ts)
    bunches.append(particle.Bunch(nPar,ts,gammas,qPb,dim))
    #print(len(bunches[-1].Pars))
for i in range(nBunch):
    for j in range(nPar):
        ts_beam[i*nPar+j] = np.random.normal(tCentroids[i],sigt,1)
        gs_beam[i*nPar+j] = np.random.normal(Gamma0,siggamma,1)
beam1 = particle.Beam(nPar,nBunch,ts_beam,gs_beam,qPb,dim)
#==============================================================================
t_start = time.time()
debug0 = 0
if debug0:
    N = int(1000)
    Vgrec = []
    Vgref = []
    for i in range(N):
        
        dt = i*tDriven0
        #cav1.PrintAll()
        for icav in range(nStation):
            cavs[i].update_Vg(dt)
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
    nRamp = int(200)
    tRec = np.zeros(nRamp)
    gRec = np.zeros(nRamp)
    tBunch = np.zeros(nBunch)
    
    dynamic_on = 0     
    rad_on = 0
    t = t0
    for i in range(nRamp):
        #print("i = ",i)
        for icav in range(nStation):
            cavs[icav].update_Vg(t)       
        for icav in range(nStation):
            cavs[icav].kick_par_beamC(beam1, dynamic_on)
            OneTurn.rad_map_par(vRad/nStation, Gamma0, beam1.gammas, rad_on)
            OneTurn.oneTurnMap_par(GMTSQ, Gamma0, beam1.gammas, beam1.ts, T0/nStation)
        t = t+T0
        tRec[i],gRec[i] = bunches[0].get_M1()
        tRec[i] -= t-t0+tCentroids[0]
        gRec[i] -=Gamma0
        for icav in range(nStation):
            cavs[icav].update_Ig(Ig,t)
    fig,axis = plt.subplots(1,1)
    axis.plot(tRec[:],gRec[:],'r.',ms=1)
    
    dynamic_on = 1
    rad_on = 1
    nTurn = int(1000)
    tRec1 = np.zeros(nTurn)
    gRec1 = np.zeros(nTurn)
    if dynamic_on and rad_on:
        for i in range(nTurn):
             if i == -100:
                 for j in range(7):
                     for k in range(nPar):
                         beam1.qs[((j)*nPar*1+k)] = 0
             #print("i = ",i)
             for icav in range(nStation):
                 cavs[icav].update_Vg(t)
             for icav in range(nStation):
                 cavs[icav].kick_par_beamC(beam1, dynamic_on)
                 OneTurn.rad_map_par(vRad/nStation, Gamma0, beam1.gammas, rad_on)
                 OneTurn.oneTurnMap_par(GMTSQ, Gamma0, beam1.gammas, beam1.ts, T0/nStation)
             t = t+T0
             tRec1[i],gRec1[i] = beam1.get_M1(0,nPar)
             tRec1[i] -= t-t0+tCentroids[0]
             for icav in range(nStation):
                 cavs[icav].update_Ig(Ig,t)
    gRec1-=Gamma0
    fig,axis = plt.subplots(1,1)
    axis.plot(tRec1[:],gRec1[:],'r.-',ms=1)
    axis.plot(tRec1[:10],gRec1[:10],'g.-',ms=5)
    axis.plot(tRec1[-10:],gRec1[-10:],'b.-',ms=5)
    
    for i in range(nBunch):
        tBunch[i] = beam1.get_M1(i,nPar)[0]-t-tCentroids[i]+t0
        #tBunch[i] = bunches[i].get_M1()[0]-t-tCentroids[i]+t0
    phiBunch = tBunch/Trf*360
    
    fig,axis = plt.subplots(1,1)
    fig.set_figheight(16)
    fig.set_figwidth(30)
    axis.plot(phiBunch,'r.',ms=10)
    axis.set_xlabel("Bunch number",fontsize=30)
    axis.set_ylabel("Phase (degree)",fontsize=30)
    axis.set_title("Phase shift along the train in first RF station.",fontsize=30)
    axis.tick_params(labelsize=30)
    
    tBunch1 = tBunch[:pattern[0]]
    tBunch2 = tBunch[pattern[2]:]
    tdiff = tBunch1-tBunch2
    fig,axis = plt.subplots(1,1)
    fig.set_figheight(16)
    fig.set_figwidth(30)
    axis.plot(tdiff[7:]*utl.c_light*1e3,'r.',ms=10)
    axis.set_xlabel("Bunch number",fontsize=30)
    axis.set_ylabel("Position Diff. [mm]", fontsize=30)
    axis.set_title("Position difference betwee two trains [mm].",fontsize=30)
    axis.tick_params(labelsize=30)
    
    save_time = time.time()
    fig.savefig("Phase_Shift_"+str(save_time)+".png",bbox_inches='tight')
    k = cav1.wrf*cav1.RoQ/2
    dtheta = 2*k*IbDC*pattern[1]*fillStep*Trf/np.abs(Vc)/np.sin(np.angle(Vc))/np.pi*180
    print("dtheta: ",dtheta)
    print("Total time : ", time.time()-t_start)
    print("Time/step : ", (time.time()-t_start)/(nTurn+nRamp))
    