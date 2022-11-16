#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 09:15:18 2022

@author: txin
"""

from setuptools import Extension, setup
from Cython.Build import cythonize
#from distutils.core import setup, Extension

import numpy

ext = [
    Extension(
        "calculate_Vg",
        ["calculate_Vg.pyx"],
        extra_compile_args=["-fopenmp"],
        extra_link_args=["-fopenmp"],
        )
    ]
setup(
    name='calculate_Vg',
    ext_modules=cythonize(ext,annotate=True),
    include_dirs=[numpy.get_include()],
    zip_safe=False,
)
