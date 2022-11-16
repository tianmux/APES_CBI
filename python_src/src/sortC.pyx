#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  6 16:10:28 2022

@author: txin
"""

cimport numpy as cnp
import numpy as np

from Cython.Compiler import Options

cdef extern from "math.h":
    double sin(double) nogil
DTYPE = np.cdouble
Options.boundscheck = False  # Deactivate bounds checking
Options.wraparound = False   # Deactivate negative indexing.
def sortC(double[:] t, double[:] gamma):
    cdef Py_ssize_t i
    cdef Py_ssize_t nPar = len(t)
    cdef double[:] tmpt = t
    cdef double[:] tmpg = gamma
    cdef double[:] 
    return t
