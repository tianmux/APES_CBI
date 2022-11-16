#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 25 11:34:51 2022

@author: txin
"""
import numpy as np
import matplotlib.pyplot as plt
import time
pi = np.pi
def Z(R, QL, wc, w,gp):
    return R/(1+gp)*w/(w+1j*QL/(1+gp)*(wc-w**2/wc))

def ImOmega_linear(M,N,r0,eta,gamma,T0,ws,w0,mu,nP,R,QL,wc,gp):
    A = M*N*r0*eta/(2*gamma*T0**2*ws)
    pw = np.array([(p*M+mu)*w0+ws for p in range(-nP,nP+1,1)])
    Zw = pw*np.real(Z(R,QL,wc,pw,gp))
    ImO = A*np.sum(Zw)

    return ImO
def plot_Z(R,QL,wc,w,gp,vline):
    Z_w = Z(R,QL,wc,w,gp)
    fig,axis = plt.subplots(1,1)
    fig.set_figheight(16)
    fig.set_figwidth(30)
    axis.plot(w/2/pi,Z_w,'r.',ms=10)
    axis.vlines(vline,0,np.max(Z_w))
    axis.vlines(vline[int(len(vline)/2)],0,np.max(Z_w),'r',linewidth=10)
    axis.set_xlabel("f",fontsize=30)
    axis.set_ylabel("Z",fontsize=30)
    axis.set_title("Z vs f",fontsize=30)

    axis.tick_params(labelsize=30)

    save_time = time.time()
    fig.savefig("Z_vs_w"+str(save_time)+".png",bbox_inches='tight')

M = 13552
qPb = 5e-10
N = qPb/1.6e-19
r0 = (1*1.6e-19)**2/(1*9.10938356e-31*3e8*3e8)
eta = 1.43e-5
gamma= 89041
T0 = 0.00033356409519815205
f0 = 1/T0
h = 216832
vSync = 0.037e9
vQuad = 0.11415341e9
V0 = (vSync**2+vQuad**2)**0.5
Ek = 45.5e9
Qs = np.sqrt(h*np.abs(V0)*eta/(2*np.pi*Ek))
ws = Qs*f0*2*np.pi
w0 = 2*np.pi*f0
mu = -2
RoQ = 6390
QL =  1e6
R = RoQ*QL
df = -2.5e3
fc = f0*h+df
wc = 2*pi*fc
gp = 0
nP = 10000
ImO = ImOmega_linear(M, N, r0, eta, gamma, T0, ws, w0, mu, nP, R, QL, wc, gp)
print(ImO)
f_step = 1
w_step = 2*pi*f_step
N_step = 20000
N_mode = 10
f_center = h*f0

w = np.array([f_center-N_step*f_step+i*f_step for i in range(2*N_step)])*2*pi
v_line = np.array([f_center+i*f0 for i in range(-N_mode,N_mode+1,1)])
plot_Z(R, QL, wc, w, gp,v_line)