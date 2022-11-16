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

#fc = 1e9
#RoQ = 50
#QL = 5e7

#df = -10
#frf = fc+df

#Ig0 = (-0.13036326764030534+0.6480610245805912j)
t0 = 0
#Ig = Ig0
#cav1 = RF.Cav(fc,frf,RoQ,QL,0)

inputs = utl.Get_Inputs('input.txt')
utl.Init_Beam(inputs)
t0s = utl.Init_Beam(inputs)
q = np.array(inputs[inputs["Name"]=="QpB"]["Values"].iloc[0].split()).astype("float")[0]
sigt = np.array(inputs[inputs["Name"]=="sigt"]["Values"].iloc[0].split()).astype("float")[0]
Gamma0 = np.array(inputs[inputs["Name"]=="Gamma0"]["Values"].iloc[0].split()).astype("float")[0]
siggamma = np.array(inputs[inputs["Name"]=="siggamma"]["Values"].iloc[0].split()).astype("float")[0]
nPar = np.array(inputs[inputs["Name"]=="Npar"]["Values"].iloc[0].split()).astype("int")[0]
qPb = np.array(inputs[inputs["Name"]=="QpB"]["Values"].iloc[0].split()).astype("float")[0]
nBunch = np.array(inputs[inputs["Name"]=="Pattern"]["Values"].iloc[0].split()).astype("int")[0]
fillStep = np.array(inputs[inputs["Name"]=="fill_step"]["Values"].iloc[0].split()).astype("int")[0]
rRing = np.array(inputs[inputs["Name"]=="R"]["Values"].iloc[0].split()).astype("float")[0]
h = np.array(inputs[inputs["Name"]=="h"]["Values"].iloc[0].split()).astype("int")[0]
RoQ = np.array(inputs[inputs["Name"]=="RoQ"]["Values"].iloc[0].split()).astype("float")[0]
QL = np.array(inputs[inputs["Name"]=="QL"]["Values"].iloc[0].split()).astype("float")[0] 
R = RoQ*QL
vRad = np.array(inputs[inputs["Name"]=="vRad"]["Values"].iloc[0].split()).astype("float")[0] 
vSync = vRad
vQuad =  np.array(inputs[inputs["Name"]=="vQuad"]["Values"].iloc[0].split()).astype("float")[0] 

T0 = 2*np.pi*rRing/utl.c_light
frf = 1/T0*h
df = -10
fc = frf+df
Trf = 1/frf
tDriven0 = 1*1/frf*h

cav1 = RF.Cav(fc,frf,RoQ,QL,0)
fL = cav1.wL/2/np.pi
Z = R/(1+1j*QL*(frf/fc-fc/frf))
Vc = vSync+1j*vQuad
IbDC = -qPb*nBunch/T0
Ig = Vc/Z-IbDC*2

cav1.Igm1 = Ig
cav1.Vgm1 = Ig*Z
cav1.Ugm1 = Ig*Z/1j/cav1.wrf
cav1.update_Ig(Ig,t0)

n_wait = int(1e9)
t0 = (n_wait)*Trf

Ib = qPb*nBunch/T0

tCentroids = np.array([t0+Trf*i*fillStep for i in range(nBunch)])


gammas = np.array(np.random.normal(Gamma0,siggamma,nPar))
dim = "2D"
#bunch1 = particle.Bunch(nPar,ts, gammas, qPb,dim)
bunches = []

for i in range(nBunch):
    ts = np.array(np.random.normal(tCentroids[i],sigt,nPar))
    ts = np.array([tCentroids[i]+0*Trf/1e2 for _ in range(nPar)])
    ts = np.sort(ts)
    bunches.append(particle.Bunch(nPar,ts,gammas,qPb,dim))
    #print(len(bunches[-1].Pars))
#==============================================================================
t_start = time.time()
debug0 = 1
if debug0:
    N = int(1e4)
    Vgrec = []
    Vgref = []
    for i in range(N):
        
        dt = i*tDriven0
        #cav1.PrintAll()
        cav1.update_Vg(dt)
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

#==============================================================================
debug1 = 0
if debug1 :
    Q = []
    Z = cav1.R/(1+1j*cav1.QL*(cav1.wrf/cav1.wL-cav1.wL/cav1.wrf))
    Ib = -q*cav1.wrf/(2*np.pi)
    Vbrec = []
    Vbref = []
    Vbref2 = []
    Vadd = -cav1.RoQ*cav1.wc/2*q
    
    for i in range(len(t0s[0])):
        if i ==0:
            Q.append(0)
        else:
            Q.append(q)
    for i in range(len(t0s[0])):
        t0 = i*tDriven0
        t = (i+1)*tDriven0
        #cav1.update_Vb(t0, t, Q,i)
        Vbrec.append(cav1.Vb+Vadd)
        Vbref.append(Ib*Z)
        Vbref2.append(Vadd/(1-np.exp((-cav1.alpha+1j*cav1.wL)*tDriven0)))
    print("Vadd : ",Vadd)
    print("VbL0: ", Vbrec[0])
    print("VbL: ", Vbrec[-1])
    print("Vbref : ",Vbref[-1])
    print("Vbref2 : ",Vbref2[-1])
    
    print("Vb dif. : ", Vbrec[-1]-Vbref[-1])
    print("Vb dif.2 : ", Vbrec[-1]-Vbref2[-1])
    
    fig,axis = plt.subplots(1,1)
    axis.plot(np.real(Vbrec[:]),'r.-')
    axis.plot(np.real(Vbref[:]),'b-')
    axis.plot(np.real(Vbref2[:]),'g-')

    # Kick to the bunch related
    #==============================================================================
    # Vg of the particles in one bunch
    
    
    cav1.kick_par(bunch1)
    Vgref_bunch = np.cos(cav1.wrf*ts-Phi)*np.abs(Ig0*Z)
    fig,axis = plt.subplots(2,1)
    axis[0].plot(ts,bunch1.Vgs,'r.')
    axis[0].plot(ts,Vgref_bunch,'b.')
    axis[1].plot(ts,Vgref_bunch-bunch1.Vgs,'g.')
    
    print("Bunch1 Vg : ",bunch1.Vgs[0])
    print("Vref : ", Vgref[-1])
    
    #==============================================================================
    # Vb of the particles in one bunch
    for i in range(len(bunch1.Pars)):
        bunch1.Vbs[i] = cav1.update_Vb(ts[i],qPb/nPar)
    fig,axis = plt.subplots(1,1)
    axis.plot(ts,np.real(bunch1.Vbs),'r.')
    
    tsamp = np.array([t0+i*Trf/1e1 for i in range(int(1e1))])
    Vbsamp = []
    Vbsamp = cav1.calcul_Vb(tsamp)
    
    fig,axis = plt.subplots(1,1)
    axis.plot(tsamp[:],np.real(Vbsamp[:]),'r.')
    #==============================================================================
    # Kick to the particle
    Vq = -cav1.RoQ*cav1.wL/2*qPb/nPar
    
    for i in range(len(bunch1.Pars)):
        bunch1.Pars[i].gamma+= (np.real(bunch1.Vgs[i]+bunch1.Vbs[i])-Vq/2)/utl.E0e
    for i in range(len(bunch1.Pars)):
        gammas[i] = bunch1.Pars[i].gamma
    
    fig,axis = plt.subplots(1,1)
    axis.plot(ts,gammas-Gamma0,'r.')

# Multi-turn kick
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
        for j in range(nBunch):
            cav1.kick_par(bunches[j],dynamic_on)
            OneTurn.rad_map(vRad,Gamma0,bunches[j],rad_on)
            OneTurn.oneTurnMap(inputs, bunches[j], T0)
        t = t+T0
        tRec[i],gRec[i] = bunches[0].get_M1()
        tRec[i] -= t-t0+tCentroids[0]
        gRec[i] -=Gamma0
        cav1.update_Ig(Ig,t)
    fig,axis = plt.subplots(1,1)
    axis.plot(tRec[:],gRec[:],'r.',ms=1)
    
    dynamic_on = 1     
    rad_on = 1
    nTurn = int(200)
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
                 OneTurn.oneTurnMap(inputs, bunches[j], T0)
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
    axis.plot(phiBunch,'r.',ms=1)
    print("Total time : ", time.time()-t_start)
    print("Time/step : ", (time.time()-t_start)/nTurn)
    