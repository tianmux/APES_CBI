#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  6 21:15:50 2022

@author: txin
"""


cimport numpy as cnp
import numpy as np

from Cython.Compiler import Options
from cython.parallel import prange

#DTYPE = np.cdouble
Options.boundscheck = False  # Deactivate bounds checking
Options.wraparound = False   # Deactivate negative indexing.
cdef extern from "complex.h":
    double complex cexp(double complex) nogil

cdef extern from "math.h":
    double sin(double) nogil
cdef extern from "math.h":
    double cos(double) nogil

def calcul_VgC(double alpha, double wrf, double wL,  double C, double L, double Tgn, double complex Ig, \
              double complex Vgm1, double complex Ugm1, double complex[::1] Vg, double[::1] t,int nThreads):
    # Solution to the RLC circuit by Laplace transformation
    cdef Py_ssize_t i
    cdef Py_ssize_t nPar = len(t)
    cdef double complex[::1] tmpVg = Vg
    cdef double[::1] tmpt = t
    
    
    cdef double complex A = alpha+1j*(wrf+wL)
    cdef double complex B = alpha+1j*(wrf-wL)
    
    for i in prange(nPar,nogil=True,num_threads=nThreads):
        #dt = tmpt[i]-Tgn
        tmpVg[i] = Ig/C*\
        (1j*wrf/A/B*cexp(1j*wrf*(tmpt[i]-Tgn))\
        +(alpha+1j*wL)/(-2j*wL*A)*cexp(-(1j*wL+alpha)*(tmpt[i]-Tgn))\
        +(alpha-1j*wL)/(2j*wL*B)*cexp((1j*wL-alpha)*(tmpt[i]-Tgn)))\
        +(Vgm1*(cos(wL*(tmpt[i]-Tgn))-alpha/wL*sin(wL*(tmpt[i]-Tgn)))\
        -Ugm1/C/L/(1*wL)*sin(wL*(tmpt[i]-Tgn)))*cexp(-alpha*(tmpt[i]-Tgn))
    
    return Vg

def calcul_Vg(double[:] alpha, double[:] wrf, double[:] wL,  double[:] C, double[:] L, double Tgn, double complex[:] Ig, \
              double complex[:] Vgm1, double complex[:] Ugm1, double complex[:,::1] Vg, double[::1] t,int nThreads):
    cdef Py_ssize_t i
    cdef Py_ssize_t nh = len(alpha)
    cdef double complex[:,::1] tmpVg = Vg
    cdef double[::1] tmpt = t
    
    for i in range(nh):
        calcul_VgC(alpha[i], wrf[i], wL[i],  C[i], L[i], Tgn, Ig[i], \
                      Vgm1[i], Ugm1[i], Vg[i], t,nThreads)
    