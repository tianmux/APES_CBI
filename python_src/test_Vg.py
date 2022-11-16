#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 15:19:57 2022

@author: txin
"""

import numpy as np

import particle
import utl
import post
import RF
import OneTurn
import matplotlib.pyplot as plt
import time
from scipy.stats import truncnorm 
import copy

#==============================================================================
t_start = time.time()

t0 = 0

inputs = utl.Get_Inputs('./inputs/input_unitTest_CEPC_Z_ZhuHang_DirectFB_CombFB.txt')
#utl.Init_Beam(inputs)
#t0s = utl.Init_Beam(inputs)
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
tau = np.array(inputs[inputs["Name"]=="tau"]["Values"].iloc[0].split()).astype("float")[0]
nWatch = np.array(inputs[inputs["Name"]=="nWatch"]["Values"].iloc[0].split()).astype("int")[0] 
sampRate = np.array(inputs[inputs["Name"]=="sampRate"]["Values"].iloc[0].split()).astype("int")[0] 

nRamp = np.array(inputs[inputs["Name"]=="nRamp"]["Values"].iloc[0].split()).astype("int")[0] 
nTrack = np.array(inputs[inputs["Name"]=="nTrack"]["Values"].iloc[0].split()).astype("int")[0] 
nSamp = int(nTrack/sampRate)+1
nThreads = np.array(inputs[inputs["Name"]=="nThreads"]["Values"].iloc[0].split()).astype("int")[0] 

# RF part
#==============================================================================
nStation =  np.array(inputs[inputs["Name"]=="nStation"]["Values"].iloc[0].split()).astype("int")[0]
h = []
RoQ = []
QL = []
vRad = []
vQuad = []
vSync = []
df = []
act = []
for i in range(nStation):
    h.append(np.array(inputs[inputs["Name"]=="h_"+str(i)]["Values"].iloc[0].split()).astype("float"))
    RoQ.append(np.array(inputs[inputs["Name"]=="RoQ_"+str(i)]["Values"].iloc[0].split()).astype("float"))
    QL.append(np.array(inputs[inputs["Name"]=="QL_"+str(i)]["Values"].iloc[0].split()).astype("float"))
    vRad.append(np.array(inputs[inputs["Name"]=="vRad_"+str(i)]["Values"].iloc[0].split()).astype("float"))
    vQuad.append(np.array(inputs[inputs["Name"]=="vQuad_"+str(i)]["Values"].iloc[0].split()).astype("float"))
    vSync.append(np.array(inputs[inputs["Name"]=="vSync_"+str(i)]["Values"].iloc[0].split()).astype("float"))
    df.append(np.array(inputs[inputs["Name"]=="df_"+str(i)]["Values"].iloc[0].split()).astype("float"))
    act.append(np.array(inputs[inputs["Name"]=="act_"+str(i)]["Values"].iloc[0].split()).astype("float"))

feed_step =  np.array(inputs[inputs["Name"]=="feed_step"]["Values"].iloc[0].split()).astype("int")[0]
gp =  np.array(inputs[inputs["Name"]=="gp"]["Values"].iloc[0].split()).astype("double")[0]
gc =  np.array(inputs[inputs["Name"]=="gc"]["Values"].iloc[0].split()).astype("double")[0]

h = np.array(h)
RoQ = np.array(RoQ)
QL = np.array(QL)
vRad = np.array(vRad)
vQuad = np.array(vQuad)
vSync = np.array(vSync)
df = np.array(df)
R = RoQ*QL
act = np.array(act)

T0 = 2*np.pi*rRing/utl.c_light
frf = 1/T0*h
damp_coeff = 2*T0/tau

fc = frf+df
Trf = 1/frf
Trf_main = Trf[0][0]
tDriven0 = 1*1/frf*h

Z = R/(1+1j*QL*(frf/fc-fc/frf))
Vc = vSync+1j*vQuad
IbDC = -qPb*nBunch/T0
Ig = (Vc/Z-IbDC*2)*act

n_wait = int(192000)
t0 = n_wait*Trf_main

N_samp = 300 # RF cycle per stride, mimic the feedback
N_step = 100 # number of strides
Nrf = int(N_samp*N_step) # Range of RF samples covered by this run
N_par = 1 # Number of samples mimic the particle between the stride
N_tot = int(N_step)*(N_par) # 

cavs = []
for i in range(nStation):
    cavs.append(RF.Cav(h[i],fc[i],frf[i],RoQ[i],QL[i],t0,N_tot,feed_step, Vc[0], Ig[0],gp,gc))


for i in range(nStation):
    
    cavs[i].Igm1 = Ig[i]
    cavs[i].Vgm1 = Ig[i]*Z[i]
    cavs[i].Ugm1 = Ig[i]*Z[i]/1j/cavs[i].wrf
    cavs[i].update_Ig(Ig[i],t0)
    #cavs[i].update_Ig(Ig[i],t0)

    cavs[i].Vg = Ig[i]*Z[i]
    
    cavs[i].update_Vg(t0+t0)
    t0 = t0+t0
    cavs[i].update_Ig(Ig[i],t0)

print(cavs[0].Vg)
print(Ig*Z)
print(cavs[0].Ug)

# Beam 
#==============================================================================
Ib = qPb*nBunch/T0

tCentroids = np.zeros(nBunch)
ibunch = 0
ibucket = 0
for i in range(nTrain):
    for j in range(pattern[i*2]):
        tCentroids[ibunch] = t0+Trf_main*ibucket 
        ibunch += 1
        ibucket += fillStep
    ibucket += fillStep*pattern[i*2+1]

gammas = np.array(np.random.normal(Gamma0,siggamma,nPar))
dim = "2D"
#bunch1 = particle.Bunch(nPar,ts, gammas, qPb,dim)
bunches = []
ts_beam = np.zeros(nBunch*nPar)
gs_beam = np.zeros(nBunch*nPar)
#for i in range(nBunch):
#    ts = truncnorm.rvs(-3, 3, size = nPar)*sigt+tCentroids[i]
#    ts = np.sort(ts)
#    bunches.append(particle.Bunch(nPar,ts,gammas,qPb,dim))
#for i in range(nBunch):
#    for j in range(nPar):
#        ts_beam[i*nPar+j] = truncnorm.rvs(-3, 3, size = 1)*sigt+tCentroids[i]
#       gs_beam[i*nPar+j] = truncnorm.rvs(-3, 3, size = 1)*siggamma+Gamma0
for i in range(nBunch):
    ts_beam[i*nPar:i*nPar+nPar] = copy.deepcopy(truncnorm.rvs(-2, 2, size = nPar)*sigt+tCentroids[i])
    gs_beam[i*nPar:i*nPar+nPar] = copy.deepcopy(truncnorm.rvs(-2, 2, size = nPar)*siggamma+Gamma0)
beam1 = particle.Beam(nPar,nBunch,ts_beam,gs_beam,qPb,dim,h[0], pattern, fillStep,nSamp)
#==============================================================================
t_start = time.time()

Trf = 1/(cavs[0].wrf[0]/np.pi/2)
dt = Trf*N_samp
ts = np.array([0.0 for i in range(N_tot)])
phis = np.array([0.0 for i in range(N_tot)])
for i in range(N_tot):
    if i%N_par == 0:
        ts[i] = float(i/N_par*dt)
    else:
        ts[i] = (int(i/N_par))*dt+Trf+Trf*(i%N_par-N_par)/N_par
        #print((i%N_par-0.5*N_par)/N_par)
for i in range(N_tot):
    phis[i] = ts[i]/Trf*2*np.pi
    
ts = ts+t0*1e9
#plt.plot(ts-t0,'r.')


Vgs = np.ndarray(N_tot).astype('complex')
tmp1 = np.ndarray(N_tot).astype('complex')
tmp2 = np.ndarray(N_tot).astype('complex')


beam1.ts = ts
cavs[0].calculate_VgC_feed(beam1, 0, N_tot,nThreads)
Vgs_cal = cavs[0].Vgs[0][:len(beam1.ts)]

#print("Tgn: ", cavs[0].Tgn)
#print("t0: ",t0)
#print("ts: ",ts)
#print(cavs[0].Vg)
#print(Ig*Z)
#print(cavs[0].Ug)

#for i in range(N_step):
    

for i in range(N_tot):
    tmp1[i],tmp2[i] = cavs[0].update_Vg(ts[i])
    #cavs[0].update_Ig(Ig*np.exp(1j*cavs[0].wrf[0]*ts[i]),ts[i])
    cavs[0].update_Ig(Ig*np.exp(1j*phis[i]),ts[i])
    #print(cavs[0].Vg[0])
    Vgs[i] = cavs[0].Vg[0]

Vgs_ref = (Ig*Z*np.exp(1j*cavs[0].wrf[0]*ts))[0]
Vgs_ref1 = (Ig*Z*(np.cos(cavs[0].wrf[0]*ts)+1j*np.sin(cavs[0].wrf[0]*ts)))[0]


fig,axis = plt.subplots(6,1)
fig.set_figheight(10)
fig.set_figwidth(20)
rng1 = 0
rng2 = rng1+100
axis[0].plot(Vgs[rng1:rng2],'r.')
axis[0].plot(Vgs_ref[rng1:rng2],'bx-')
axis[0].legend(["Vg_from_update_Vg","Vg_from_pure_exp"])

axis[1].plot(Vgs_cal[rng1:rng2],'r.')
axis[1].plot(Vgs_ref[rng1:rng2],'bx-')
axis[1].legend(["Vg_from_cal_Vg","Vg_from_pure_exp"])

axis[2].plot(Vgs_ref[rng1:rng2],'r.')
axis[2].plot(Vgs_ref1[rng1:rng2],'bx-')
axis[2].legend(["Vg_from_pure_exp","Vg_from_pure_cos"])

axis[3].plot(((Vgs_ref-Vgs_ref1)/Vgs_ref)[rng1:rng2],'b-')
axis[3].legend(["Vgs_ref-Vgs_ref1)/Vgs_ref"])

axis[4].plot(((Vgs_ref-Vgs)/np.abs(Vgs_ref))[rng1:rng2],'b-')
axis[4].legend(["(Vgs_ref-Vgs)/Vgs_ref"])

axis[5].plot(((Vgs_ref-Vgs_cal)/np.abs(Vgs_ref))[rng1:rng2],'b-')
axis[5].legend(["(Vgs_ref-Vgs_cal)/Vgs_ref"])

fig,axis = plt.subplots(4,1)
fig.set_figheight(10)
fig.set_figwidth(10)
rng1 = 0
rng2 = rng1+100
axis[0].plot(tmp1[rng1:rng2],'r.-')
axis[0].legend(["Ig part"])

axis[1].plot(tmp2[rng1:rng2],'r.-')
axis[1].legend(["Vg and Ug part"])

axis[2].plot((tmp1+tmp2)[rng1:rng2],'r.-')
axis[2].legend(["tmp1+tmp2"])

axis[3].plot((tmp1+tmp2-Vgs_ref)[rng1:rng2],'r.-')
axis[3].legend(["tmp1+tmp2-Vgs_ref"])

