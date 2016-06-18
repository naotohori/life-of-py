#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2016/06/18
@author: Naoto Hori
'''

import sys
import math
import numpy as np
from cafysis.lib_f2py import py_distance2_hist_nt
from cafysis.file_io.dcd import DcdFile
from cafysis.file_io.psf import PsfFile

if len(sys.argv) != 4:
    print ('\n Usage: SCRIPT [input DCD] [input PSF] [output]\n')
    sys.exit(2)
    
psffile = PsfFile(sys.argv[2])
psffile.open_to_read()
psf = psffile.read_all()

dcd = DcdFile(sys.argv[1])

id0_Mg = []
id0_P  = []
id0_K  = []
id0_Cl = []
for i,a in enumerate(psf.atoms):
    if a.atom_name.strip() == 'P':
        id0_P.append(i)
    elif a.atom_name.strip() == 'Mg':
        id0_Mg.append(i)
    elif a.atom_name.strip() == 'K':
        id0_K.append(i)
    elif a.atom_name.strip() == 'Cl':
        id0_Cl.append(i)
nP = len(id0_P)

#print id0_Mg
#print id0_P

f_out = open(sys.argv[-1],'w')

r_hist = 20.0
r2_hist = r_hist * r_hist
dr = 0.1
bins = [(x*dr)**2 for x in range(200+1)]
nbin = len(bins) - 1

hist_Mg = np.zeros((nP,nbin))
hist_K  = np.zeros((nP,nbin))
hist_Cl = np.zeros((nP,nbin))

dcd.open_to_read()
dcd.read_header()

nstep = 0
while dcd.has_more_data() :
    nstep += 1

    data = dcd.read_onestep_npF()

    if len(id0_Mg) > 0:
        dist2, idx = py_distance2_hist_nt.distance2_hist_nt( data, id0_P, id0_Mg, r2_hist)
        
        for iP in range(nP):
            hist, B = np.histogram(dist2[:idx[iP],iP], bins=bins)
            hist_Mg[iP,:] += hist[:]

    if len(id0_K) > 0:
        dist2, idx = py_distance2_hist_nt.distance2_hist_nt( data, id0_P, id0_K, r2_hist)
        
        for iP in range(nP):
            hist, B = np.histogram(dist2[:idx[iP],iP], bins=bins)
            hist_K[iP,:] += hist[:]

    if len(id0_Cl) > 0:
        dist2, idx = py_distance2_hist_nt.distance2_hist_nt( data, id0_P, id0_Cl, r2_hist)

        for iP in range(nP):
            hist, B = np.histogram(dist2[:idx[iP],iP], bins=bins)
            hist_Cl[iP,:] += hist[:]

dcd.close()

for iP in range(nP):

    f_out.write('#%5i\n' % (iP+1))

    if len(id0_Mg) > 0 and len(id0_K) > 0:
        for i in range(nbin):
            r1 = math.sqrt(bins[i])
            r2 = math.sqrt(bins[i+1])
            div_factor = float(nstep) * 4.0 / 3.0 * math.pi * (r2**3 - r1**3)
            f_out.write('%4.1f %4.1f %7.5f %7.5f %7.5f\n' % (r1,r2, 
                                         hist_Mg[iP,i] / div_factor, 
                                         hist_K[iP,i]  / div_factor,
                                         hist_Cl[iP,i] / div_factor) )
    elif len(id0_K) > 0:
        for i in range(nbin):
            r1 = math.sqrt(bins[i])
            r2 = math.sqrt(bins[i+1])
            div_factor = float(nstep) * 4.0 / 3.0 * math.pi * (r2**3 - r1**3)
            f_out.write('%4.1f %4.1f %7.5f %7.5f %7.5f\n' % (r1,r2,
                                         0.0, 
                                         hist_K[iP,i] / div_factor, 
                                         hist_Cl[iP,i]/ div_factor ) )
    f_out.write('\n\n')

f_out.close()

'''
! dist2_m,idx_m = distance2_hist_nt(xyz,id0_p,id0_m,r2_hist,[nmp,np,nm])
! 
! Wrapper for ``distance2_hist_nt``.
! 
! Parameters
! ----------
! xyz : input rank-2 array('d') with bounds (3,nmp)
! id0_p : input rank-1 array('i') with bounds (np)
! id0_m : input rank-1 array('i') with bounds (nm)
! r2_hist : input float
! 
! Other Parameters
! ----------------
! nmp : input int, optional
!     Default: shape(xyz,1)
! np : input int, optional
!     Default: len(id0_p)
! nm : input int, optional
!     Default: len(id0_m)
! 
! Returns
! -------
! dist2_m : rank-2 array('d') with bounds (nm,np)
! idx_m : rank-1 array('i') with bounds (np)
'''
