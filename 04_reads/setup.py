#!/bin/python3
#-*- coding:utf-8 -*-
#cython: language_level=3

from distutils.core import Extension, setup
from Cython.Build import cythonize

names=["config", "ahocorasick", "read", "reads"]

for name in names :
    ext = Extension(name=name, sources=[f"{name}.pyx"])
    setup(ext_modules=cythonize(ext), script_args=["build_ext", "--inplace"])