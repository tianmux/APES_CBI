#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 11:59:36 2022

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

# ==============================================================================
t_start = time.time()

t0 = 0

inputs = utl.Get_Inputs('./inputs/input_unitTest_CEPC_Higgs_DFB_OTFB.txt')
# utl.Init_Beam(inputs)
#t0s = utl.Init_Beam(inputs)
q = np.array(inputs[inputs["Name"] == "QpB"]
             ["Values"].iloc[0].split()).astype("float")[0]
sigt = np.array(inputs[inputs["Name"] == "sigt"]
                ["Values"].iloc[0].split()).astype("float")[0]
Gamma0 = np.array(inputs[inputs["Name"] == "Gamma0"]
                  ["Values"].iloc[0].split()).astype("float")[0]
GMTSQ = np.array(inputs[inputs["Name"] == "GMTSQ"]
                 ["Values"].iloc[0].split()).astype("float")[0]
siggamma = np.array(inputs[inputs["Name"] == "siggamma"]
                    ["Values"].iloc[0].split()).astype("float")[0]
nPar = np.array(inputs[inputs["Name"] == "Npar"]
                ["Values"].iloc[0].split()).astype("int")[0]
qPb = np.array(inputs[inputs["Name"] == "QpB"]
               ["Values"].iloc[0].split()).astype("float")[0]
nTrain = np.array(inputs[inputs["Name"] == "nTrain"]
                  ["Values"].iloc[0].split()).astype("int")[0]
pattern = np.array(inputs[inputs["Name"] == "Pattern"]
                   ["Values"].iloc[0].split()).astype("int")

nBunch = np.array(inputs[inputs["Name"] == "nBunch"]
                  ["Values"].iloc[0].split()).astype("int")[0]
fillStep = np.array(inputs[inputs["Name"] == "fill_step"]
                    ["Values"].iloc[0].split()).astype("int")[0]

rRing = np.array(inputs[inputs["Name"] == "R"]
                 ["Values"].iloc[0].split()).astype("float")[0]
tau = np.array(inputs[inputs["Name"] == "tau"]
               ["Values"].iloc[0].split()).astype("float")[0]
nWatch = np.array(inputs[inputs["Name"] == "nWatch"]
                  ["Values"].iloc[0].split()).astype("int")[0]
sampRate = np.array(inputs[inputs["Name"] == "sampRate"]
                    ["Values"].iloc[0].split()).astype("int")[0]

nRamp = np.array(inputs[inputs["Name"] == "nRamp"]
                 ["Values"].iloc[0].split()).astype("int")[0]
nTrack = np.array(inputs[inputs["Name"] == "nTrack"]
                  ["Values"].iloc[0].split()).astype("int")[0]
nSamp = int(nTrack/sampRate)+1
nThreads = np.array(inputs[inputs["Name"] == "nThreads"]
                    ["Values"].iloc[0].split()).astype("int")[0]

# RF part
# ==============================================================================
nStation = np.array(inputs[inputs["Name"] == "nStation"]
                    ["Values"].iloc[0].split()).astype("int")[0]
h = []
RoQ = []
QL = []
vRad = []
vQuad = []
vSync = []
df = []
act = []
for i in range(nStation):
    h.append(np.array(inputs[inputs["Name"] == "h_"+str(i)]
             ["Values"].iloc[0].split()).astype("float"))
    RoQ.append(np.array(inputs[inputs["Name"] == "RoQ_"+str(i)]
               ["Values"].iloc[0].split()).astype("float"))
    QL.append(np.array(inputs[inputs["Name"] == "QL_"+str(i)]
              ["Values"].iloc[0].split()).astype("float"))
    vRad.append(np.array(inputs[inputs["Name"] == "vRad_"+str(i)]
                ["Values"].iloc[0].split()).astype("float"))
    vQuad.append(np.array(
        inputs[inputs["Name"] == "vQuad_"+str(i)]["Values"].iloc[0].split()).astype("float"))
    vSync.append(np.array(
        inputs[inputs["Name"] == "vSync_"+str(i)]["Values"].iloc[0].split()).astype("float"))
    df.append(np.array(inputs[inputs["Name"] == "df_"+str(i)]
              ["Values"].iloc[0].split()).astype("float"))
    act.append(np.array(inputs[inputs["Name"] == "act_"+str(i)]
               ["Values"].iloc[0].split()).astype("float"))

feed_step = np.array(inputs[inputs["Name"] == "feed_step"]
                     ["Values"].iloc[0].split()).astype("int")[0]
gp = np.array(inputs[inputs["Name"] == "gp"]
              ["Values"].iloc[0].split()).astype("double")[0]
gc = np.array(inputs[inputs["Name"] == "gc"]
              ["Values"].iloc[0].split()).astype("double")[0]
epsilon = np.array(inputs[inputs["Name"] == "epsilon"]
              ["Values"].iloc[0].split()).astype("double")[0]
delay_dt = np.array(inputs[inputs["Name"] == "delay_dt"]
              ["Values"].iloc[0].split()).astype("double")[0]
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
lmd = utl.c_light/frf
damp_coeff = 2*T0/tau
excite_coeff = siggamma*np.sqrt(3.0)*np.sqrt(4*T0/tau)
fc = frf+df
Trf = 1/frf
Trf_main = Trf[0][0]
tDriven0 = 1*1/frf*h

Z = R/(1+1j*QL*(frf/fc-fc/frf))
Vc = vSync+1j*vQuad
IbDC = -qPb*nBunch/T0
Ig = (Vc/Z-IbDC*2)*act
print("Igref = ",Ig)
n_wait = int(1e9)
t0 = n_wait*Trf_main

cavs = []
for i in range(nStation):
    cavs.append(RF.Cav(h[i], fc[i], frf[i], RoQ[i], QL[i],
                t0, nPar*nBunch, feed_step, Vc[0], Ig[0], gp, gc,epsilon,delay_dt))


for i in range(nStation):
    
    cavs[i].Igm1 = Ig[i]
    cavs[i].Vgm1 = Ig[i]*Z[i]
    cavs[i].Ugm1 = Ig[i]*Z[i]/1j/cavs[i].wrf
    cavs[i].update_Ig(Ig[i],t0)
    #cavs[i].update_Ig(Ig[i],t0)

    cavs[i].Vg = Ig[i]*Z[i]
    
    cavs[i].update_Vg(t0+t0)    
    cavs[i].update_feed_times(t0)

    t0 = t0+t0
    cavs[i].update_Ig(Ig[i],t0)

#n_wait = int(1e9)
#t0 = n_wait*Trf_main

# Beam
# ==============================================================================
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

gammas = np.array(np.random.normal(Gamma0, siggamma, nPar))
dim = "2D"
#bunch1 = particle.Bunch(nPar,ts, gammas, qPb,dim)
bunches = []
ts_beam = np.zeros(nBunch*nPar)
gs_beam = np.zeros(nBunch*nPar)
# for i in range(nBunch):
#    ts = truncnorm.rvs(-3, 3, size = nPar)*sigt+tCentroids[i]
#    ts = np.sort(ts)
#    bunches.append(particle.Bunch(nPar,ts,gammas,qPb,dim))
# for i in range(nBunch):
#    for j in range(nPar):
#        ts_beam[i*nPar+j] = truncnorm.rvs(-3, 3, size = 1)*sigt+tCentroids[i]
#       gs_beam[i*nPar+j] = truncnorm.rvs(-3, 3, size = 1)*siggamma+Gamma0
for i in range(nBunch):
    ts_beam[i*nPar:i*nPar +
            nPar] = copy.deepcopy(truncnorm.rvs(-2, 2, size=nPar)*sigt+tCentroids[i])
    gs_beam[i*nPar:i*nPar +
            nPar] = copy.deepcopy(truncnorm.rvs(-2, 2, size=nPar)*siggamma+Gamma0)
beam1 = particle.Beam(nPar, nBunch, ts_beam, gs_beam, qPb,
                      dim, h[0], pattern, fillStep, nSamp)
# ==============================================================================
t_start = time.time()

# Multi-turn kick

debug2 = 1
if debug2:

    tRec = np.zeros(nRamp)
    gRec = np.zeros(nRamp)
    tBunch = np.zeros(nBunch)

    dynamic_on = 0
    rad_on = 0
    feed_on = 0

    t = t0
    for i in range(nRamp):
        if nRamp >= 10 and i % int(nRamp/10) == 0:
            print("Ramping ", int(i*100/nRamp), "%.")
        for icav in range(nStation):
            cavs[icav].update_Vg(t)
        for icav in range(nStation):
            if feed_on == 1:
                cavs[icav].kick_par_beamC_feed(
                    beam1, dynamic_on, feed_on, nThreads)
            if feed_on == 0:
                cavs[icav].kick_par_beamC(beam1, dynamic_on, nThreads)

            OneTurn.rad_map_beam(
                vRad[0], Gamma0, beam1.gammas, rad_on, damp_coeff,excite_coeff)
            OneTurn.oneTurnMap_beam(
                GMTSQ, Gamma0, beam1.gammas, beam1.ts, T0/nStation, dynamic_on)
            cavs[icav].update_feed_times(T0)
        t = t+T0
        tRec[i], gRec[i] = beam1.get_M1(0, nPar)
        tRec[i] -= t-t0+tCentroids[0]
        gRec[i] -= Gamma0
        # for icav in range(nStation):
        # cavs[icav].update_Ig(Ig[icav]*np.exp(1j*cavs[icav].wrf*t),t)
    fig, axis = plt.subplots(1, 1)
    axis.plot(tRec[:], gRec[:], 'r.', ms=1)

    dynamic_on = 1
    rad_on = 1
    feed_on = 1
    tRec1 = np.zeros(nTrack)
    gRec1 = np.zeros(nTrack)

    if dynamic_on or rad_on:
        for i in range(nTrack):
            if nTrack >= 10 and i % int(nTrack/10) == 0:
                print("Tracking ", int(i*100/nTrack), "%.")
            if i % sampRate == 0:
                beam1.get_M1s(int(i/sampRate), Trf, t, t0, tCentroids)
                #print("Ig_last = ",cavs[0].Ig)


            if i == -100:
                for j in range(7):
                    for k in range(nPar):
                        beam1.qs[((j+120)*nPar*1+k)] = 0
            # for icav in range(nStation):
                # cavs[icav].update_Vg(t)
            for icav in range(nStation):
                if feed_on == 1:
                    cavs[icav].kick_par_beamC_feed(
                        beam1, dynamic_on, feed_on, nThreads)
                if feed_on == 0:
                    cavs[icav].kick_par_beamC(beam1, dynamic_on, nThreads)

                OneTurn.rad_map_beam(vRad[0], Gamma0, beam1.gammas, rad_on, damp_coeff,excite_coeff)
                OneTurn.oneTurnMap_beam(GMTSQ, Gamma0, beam1.gammas, beam1.ts, T0/nStation, dynamic_on)
                cavs[icav].update_feed_times(T0)
            t = t+T0
            tRec1[i], gRec1[i] = beam1.get_M1(0, nPar)
            tRec1[i] -= t-t0+tCentroids[0]
            # for icav in range(nStation):
            # cavs[icav].update_Ig(Ig[icav],t)
        beam1.get_M1s(nSamp-1, Trf, t, t0, tCentroids)
    fig, axis = plt.subplots(2, 1)
    axis[0].plot(cavs[0].Vgs[0][0:15].real, '.-')
    axis[1].plot(cavs[0].Vbs[0][0:15].real, '.-')

    fig, axis = plt.subplots(2, 1)
    axis[0].plot(cavs[0].Vgs[0].real, '.-')
    axis[1].plot(cavs[0].Vbs[0].real, '.-')
    print(cavs[0].Vgs[0][-1].real)
    gRec1 -= Gamma0
    fig, axis = plt.subplots(1, 1)
    axis.plot(tRec1[:120], gRec1[:120], 'r.-', ms=1)
    axis.plot(tRec1[:10], gRec1[:10], 'g.-', ms=5)
    axis.plot(tRec1[-10:], gRec1[-10:], 'b.-', ms=5)

    for i in range(nBunch):
        #print('test1.1')
        tBunch[i] = beam1.get_M1(i, nPar)[0]-t-tCentroids[i]+t0
        #print('test1.2')
    phiBunch = tBunch/Trf_main*360

    post.plot_controids(beam1, sampRate, nTrack-1,Gamma0)
    fig, axis = plt.subplots(1, 1)
    fig.set_figheight(16)
    fig.set_figwidth(30)
    axis.plot(phiBunch, 'r.', ms=10)
    axis.set_xlabel("Bunch number", fontsize=30)
    axis.set_ylabel("Phase (degree)", fontsize=30)
    axis.set_title(
        "Phase shift along the train in first RF station.", fontsize=30)
    axis.tick_params(labelsize=30)
    print("Phase difference : ",phiBunch[0]-phiBunch[-1])
    post.plot_phi_diff(phiBunch, pattern)
    post.plot_z_diff(phiBunch, pattern,lmd[0][0])

    #tBunch1 = tBunch[:pattern[0]]
    #tBunch2 = tBunch[pattern[2]:]
    #tdiff = tBunch1-tBunch2
    #fig,axis = plt.subplots(1,1)
    # fig.set_figheight(16)
    # fig.set_figwidth(30)
    # axis.plot(tdiff[7:]*utl.c_light*1e3,'r.',ms=10)
    #axis.set_xlabel("Bunch number",fontsize=30)
    #axis.set_ylabel("Position Diff. [mm]", fontsize=30)
    #axis.set_title("Position difference betwee two trains [mm].",fontsize=30)
    # axis.tick_params(labelsize=30)

    save_time = time.time()
    fig.savefig("Phase_Shift_"+str(save_time)+".png", bbox_inches='tight')
    k = cavs[0].wrf*cavs[0].RoQ/2
    dtheta = 2*k*IbDC*pattern[1]*fillStep*Trf_main / \
        np.abs(Vc[0][0])/np.sin(np.angle(Vc[0][0]))/np.pi*180
    dz = dtheta/360*Trf_main*utl.c_light*1000
    
    Amp = post.get_mode_amp(beam1, Gamma0)
    istart = 5
    iend = 30
    istep = 1
    b_gues = 1
    mus = np.array([-1, -2, -3, -4, -5, -6,-7,-8,-9])
    OmegaIms = np.zeros(len(mus))
    mu_idx = 0
    for mu in mus:
        #print("mu = ",mus[mu_idx])
        popt, pcov = post.fit_growth_rate(np.abs(Amp), mu, istart, iend, istep, b_gues)
        rng, Amp_fit = post.get_fitted_curve(np.abs(Amp), mu, istart, iend, istep, popt)
        post.plot_mode_amp(np.abs(Amp), Amp_fit, rng, mu, istart, iend)
        OmegaIms[mu_idx] = popt[1]/T0
        mu_idx += 1

    #fig, axis = plt.subplots(1, 1)
    # fig.set_figheight(16)
    # fig.set_figwidth(30)
    #axis.plot(rng, Amp.T[mu][istart:iend], 'r.', ms=10)
    #axis.plot(rng, Amp_fit, 'b-', ms=10)

    #axis.set_xlabel("Turn", fontsize=30)
    #axis.set_ylabel("Amp", fontsize=30)
    #xis.set_title("Amplitude of mode "+str(mu), fontsize=30)
    # xis.tick_params(labelsize=30)
    # plt.show()
    # print(popt[1]/T0)
    #mus = np.array([i+1 for i in range(9)])
    #OmegaIms = post.get_growth_rate(Amp, mus, istart, iend, istep, b_gues)/T0
    print(OmegaIms)

    print("dtheta: ", dtheta)
    print("dz: ", dz)
    print("Total time : ", time.time()-t_start)
    print("Time/step : ", (time.time()-t_start)/(nTrack+nRamp))
    fn = "Imag_Omega"+str(time.time())+"_gp="+str(gp)+"_gc="+str(gc)+".txt"
    post.write_output(fn, OmegaIms, mus)