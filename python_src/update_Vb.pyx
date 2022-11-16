#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 09:50:57 2022

@author: txin
"""
cimport numpy as np
import numpy as np

from Cython.Compiler import Options
from cython.parallel import prange

cdef extern from "complex.h":
    double complex cexp(double complex) nogil

DTYPE = np.cdouble
Options.boundscheck = False  # Deactivate bounds checking
Options.wraparound = False   # Deactivate negative indexing.
cdef void update_Vb(complex[:] Vb, double[:] t, double[:] Vadd, double wrf, double alpha, complex Vb0, double t0) nogil:
    cdef Py_ssize_t i
    cdef Py_ssize_t nPar = len(Vb)
    cdef double complex[:] tmpVb = Vb
    cdef double[:] tmpt = t
    cdef double[:] tmpVadd = Vadd
    cdef double complex iwrf = 1j*wrf-alpha
    
    tmpVb[0] = Vb0*cexp((1j*wrf-alpha)*(tmpt[0]-t0))+tmpVadd[0]
    for i in range(nPar-1):
        tmpVb[i+1] = tmpVb[i]*cexp((1j*wrf-alpha)*(tmpt[i+1]-tmpt[i]))+tmpVadd[i+1]
    return 

def update_Vb_beam(complex[:,::1] Vb, double[:] t, double[:,::1] Vadd, double[:] wrf, double[:] alpha, complex[:] Vb0, double t0, int nThreads):
    cdef Py_ssize_t nh = len(alpha)
    cdef Py_ssize_t i
    for i in prange(nh,nogil=True,num_threads=nThreads):
        update_Vb(Vb[i],t,Vadd[i],wrf[i],alpha[i],Vb0[i],t0)
        