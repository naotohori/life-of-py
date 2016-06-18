#!/bin/sh
#f2py -c -m py_bestfit --fcompiler=intel --compiler=intel bestfit.F90 dsvdc.f
#f2py -c -m py_bestfit --f90exec=/opt/intel/bin/ifort bestfit.F90 dsvdc.f

# for 64-bit
#f2py -c -m py_bestfit --fcompiler=intelem --compiler=intelem bestfit.F90 dsvdc.f

# Mac
f2py -c -m py_drid --fcompiler=gnu95 drid.F90
f2py -c -m py_ddrid --fcompiler=gnu95 ddrid.F90
f2py -c -m py_distance2_hist --fcompiler=gnu95 distance2_hist.F90
f2py -c -m py_count_bound --fcompiler=gnu95 count_bound.F90
f2py -c -m py_dcd_r2_histogram --fcompiler=gnu95 dcd_r2_histogram.F90
f2py -c -m py_distance2_hist_nt --fcompiler=gnu95 distance2_hist_nt.F90
