#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 09:15:18 2022

@author: txin
"""

from setuptools import setup
from Cython.Build import cythonize

import numpy

setup(
    name='Hello world app',
    ext_modules=cythonize("update_Vb.pyx",annotate=True),
    include_dirs=[numpy.get_include()],
    zip_safe=False,
)
