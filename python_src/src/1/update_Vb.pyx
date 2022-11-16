#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 09:50:57 2022

@author: txin
"""
import numpy as np

def update_Vb_beam(Vb, dt, Vadd, wrf, alpha):
    Vb[0] = Vb[-1]*np.exp((1j*wrf-alpha)*(dt[0]))+Vadd[0]
    print(Vb[0])
    for i in range(len(Vb)-1):
        Vb[i+1] = Vb[i]*np.exp((1j*wrf-alpha)*(dt[i+1]))+Vadd[i+1]
    return Vb
