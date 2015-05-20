#!/bin/sh
#f2py -c -m py_bestfit --fcompiler=intel --compiler=intel bestfit.F90 dsvdc.f
#f2py -c -m py_bestfit --f90exec=/opt/intel/bin/ifort bestfit.F90 dsvdc.f

# for 64-bit
#f2py -c -m py_bestfit --fcompiler=intelem --compiler=intelem bestfit.F90 dsvdc.f

# Mac
f2py -c -m py_drid --fcompiler=gnu95 drid.F90
