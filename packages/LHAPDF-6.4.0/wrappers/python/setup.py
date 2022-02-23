#! /usr/bin/env python

import os
from distutils.core import setup
from glob import glob
from distutils.extension import Extension

incdir_src = os.path.abspath("../../include")
incdir_build = os.path.abspath("../../include")
libdir = os.path.abspath("../../src/.libs")

## Configure the C++ extension and LHAPDF package
ext = Extension("lhapdf",
                ["lhapdf.cpp"],
                include_dirs=[incdir_src, incdir_build],
                extra_compile_args=["-I/scratch/cb27g11/THDM_T3PS_scanner/packages/lhapdf-6.4.0/include"],
                library_dirs=[libdir],
                language="C++",
                libraries=["stdc++", "LHAPDF"])
setup(name="LHAPDF",
      version="6.4.0",
      ext_modules=[ext])
