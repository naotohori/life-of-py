#!/bin/sh

CF2PYC='gnu'

#CF2PYF='gnu95'
CF2PYF='gfortran'

#f2py -c -m py_bestfit --fcompiler=intel --compiler=intel bestfit.F90 dsvdc.f
#f2py -c -m py_bestfit --f90exec=/opt/intel/bin/ifort bestfit.F90 dsvdc.f

# for 64-bit
#f2py -c -m py_bestfit --fcompiler=intelem --compiler=intelem bestfit.F90 dsvdc.f
#f2py -c -m py_drid --fcompiler=intelem --compiler=intelem drid.F90
#f2py -c -m py_ddrid --fcompiler=intelem --compiler=intelem ddrid.F90
#f2py -c -m py_distance2_hist --fcompiler=intelem --compiler=intelem distance2_hist.F90
#f2py -c -m py_count_bound --fcompiler=intelem --compiler=intelem count_bound.F90
#f2py -c -m py_count_bound_total --fcompiler=intelem --compiler=intelem count_bound_total.F90
#f2py -c -m py_dcd_r2_histogram --fcompiler=intelem --compiler=intelem dcd_r2_histogram.F90
#f2py -c -m py_distance2_hist_nt --fcompiler=intelem --compiler=intelem distance2_hist_nt.F90
#f2py -c -m py_distance_count_within --fcompiler=intelem --compiler=intelem distance_count_within.F90
#f2py -c -m py_contact --fcompiler=intelem --compiler=intelem contact.F90

# Mac
f2py -c -m py_bestfit --fcompiler=$CF2PYF bestfit.F90 dsvdc.f
f2py -c -m py_drid --fcompiler=$CF2PYF drid.F90
f2py -c -m py_ddrid --fcompiler=$CF2PYF ddrid.F90
f2py -c -m py_distance2_hist --fcompiler=$CF2PYF distance2_hist.F90
f2py -c -m py_count_bound --fcompiler=$CF2PYF count_bound.F90
f2py -c -m py_dcd_r2_histogram --fcompiler=$CF2PYF dcd_r2_histogram.F90
f2py -c -m py_distance2_hist_nt --fcompiler=$CF2PYF distance2_hist_nt.F90
f2py -c -m py_distance_count_within --fcompiler=$CF2PYF distance_count_within.F90
