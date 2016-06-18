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
bins = [(x*0.1)**2 for x in range(200+1)]
nbin = len(bins) - 1
#hist_Mg = [np.zeros( (nbin,)) * nP]
#hist_K  = [np.zeros( (nbin,)) * nP]
#hist_Cl = [np.zeros( (nbin,)) * nP]
dist2_Mg = []
dist2_K = []
dist2_Cl = []

dcd.open_to_read()
dcd.read_header()

nstep = 0
while dcd.has_more_data() :
    nstep += 1

    data = dcd.read_onestep_npF()

    if len(id0_Mg) > 0:
        dist2, idx = py_distance2_hist.distance2_hist_nt( data, id0_P, id0_Mg, r2_hist)
        
        for i in range(nP):
            dist2_Mg += dist2[iP][:idx[iP]]

    if len(id0_K) > 0:
        dist2, idx = py_distance2_hist.distance2_hist_nt( data, id0_P, id0_K, r2_hist)
        
        for i in range(nP):
            dist2_K += dist2[iP][:idx[iP]]

    if len(id0_Cl) > 0:
        dist2, idx = py_distance2_hist.distance2_hist_nt( data, id0_P, id0_Cl, r2_hist)

        for i in range(nP):
            dist2_Cl += dist2[iP][:idx[iP]]

dcd.close()

# B is dummy
hist_Mg, B = np.histogram(dist2_Mg, bins=bins)
hist_K,  B = np.histogram(dist2_K, bins=bins)
hist_Cl, B = np.histogram(dist2_Cl, bins=bins)

n_hist_Mg = len(dist2_Mg)
n_hist_K  = len(dist2_K)
n_hist_Cl = len(dist2_Cl)

if len(id0_Mg) > 0 and len(id0_K) > 0:
    for i in range(nbin):
        f_out.write('%f %f %i %f %i %f %i %f\n' % (math.sqrt(bins[i]),math.sqrt(bins[i+1]), 
                                         hist_Mg[i], hist_Mg[i] / float(n_hist_Mg),
                                         hist_K[i],  hist_K[i]  / float(n_hist_K),
                                         hist_Cl[i], hist_Cl[i] / float(n_hist_Cl),) )
elif len(id0_K) > 0:
    for i in range(nbin):
        f_out.write('%f %f %i %f %i %f %i %f\n' % (math.sqrt(bins[i]),math.sqrt(bins[i+1]), 
                                         0, 0.0,
                                         hist_K[i],  hist_K[i]  / float(n_hist_K),
                                         hist_Cl[i], hist_Cl[i] / float(n_hist_Cl),) )

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
