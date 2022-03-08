#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2015/06/03
@author: Naoto Hori
'''
import sys
import math
import numpy as np
from lop.lib_f2py import py_distance2_hist
from lop.file_io.dcd import DcdFile
from lop.file_io.psf import PsfFile

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

#print id0_Mg
#print id0_P

f_out = open(sys.argv[-1],'w')

r_hist = 20.0
r2_hist = r_hist * r_hist
bins = [(x*0.1)**2 for x in range(200+1)]
nbin = len(bins) - 1
hist_Mg = np.zeros( (nbin,))
hist_K  = np.zeros( (nbin,))
hist_Cl = np.zeros( (nbin,))
hist_P  = np.zeros( (nbin,))

dcd.open_to_read()
dcd.read_header()

nstep = 0
while dcd.has_more_data() :
    nstep += 1

    data = dcd.read_onestep_npF()

    if len(id0_Mg) > 0:
        dist2_Mg, idx_Mg, = py_distance2_hist.distance2_hist( data, id0_P, id0_Mg, r2_hist)

        H, B = np.histogram(dist2_Mg[:idx_Mg], bins=bins)
        hist_Mg = hist_Mg + H

    if len(id0_K) > 0:
        dist2_K, idx_K = py_distance2_hist.distance2_hist( data, id0_P, id0_K, r2_hist)

        H, B = np.histogram(dist2_K[:idx_K], bins=bins)
        hist_K = hist_K + H

    if len(id0_Cl) > 0:
        dist2_Cl, idx_Cl = py_distance2_hist.distance2_hist( data, id0_P, id0_Cl, r2_hist)

        H, B = np.histogram(dist2_Cl[:idx_Cl], bins=bins)
        hist_Cl = hist_Cl + H

    if len(id0_P) > 0:
        dist2_P, idx_P = py_distance2_hist.distance2_hist( data, id0_P, id0_P, r2_hist)

        H, B = np.histogram(dist2_P[:idx_P], bins=bins)
        hist_P = hist_P + H
dcd.close()

factor = float(len(id0_P)) * float(nstep) * 4.0 / 3.0 * math.pi

f_out.write('#nstep: %i\n#\n' % (nstep,))
f_out.write('#r1 r2 hist_Mg rho(r)_Mg hist_K rho(r)_K hist_Cl rho(r)_Cl hist_P rho(r)_P\n')
f_out.write('## rho(r1<r<r2) = hist / (float(#P) * float(nstep) * 4.0 / 3.0 * math.pi * (r2**3 - r1**3))\n')
if len(id0_Mg) > 0 and len(id0_K) > 0:
    for i in range(nbin):
        r1 = math.sqrt(bins[i])
        r2 = math.sqrt(bins[i+1])
        div_factor = factor * (r2**3 - r1**3)
        f_out.write('%f %f  %i %g  %i %g  %i %g  %i %g\n' % (r1,r2, 
                                         hist_Mg[i], hist_Mg[i] / div_factor,
                                         hist_K[i],  hist_K[i]  / div_factor,
                                         hist_Cl[i], hist_Cl[i] / div_factor,
                                         hist_P[i],  hist_P[i] / div_factor,) )
elif len(id0_K) > 0:
    for i in range(nbin):
        r1 = math.sqrt(bins[i])
        r2 = math.sqrt(bins[i+1])
        div_factor = factor * (r2**3 - r1**3)
        f_out.write('%f %f  %i %g  %i %g  %i %g  %i %g\n' % (r1,r2,
                                         0, 0.0,
                                         hist_K[i],  hist_K[i]  / div_factor,
                                         hist_Cl[i], hist_Cl[i] / div_factor,
                                         hist_P[i],  hist_P[i] / div_factor,) )

f_out.close()


'''
dist2_m,idx_m = distance2_hist(xyz,id0_p,id0_m,r2_hist,[nmp,np,nm])

Wrapper for ``distance2_hist``.

Parameters
----------
xyz : input rank-2 array('d') with bounds (3,nmp)
id0_p : input rank-1 array('i') with bounds (np)
id0_m : input rank-1 array('i') with bounds (nm)
r2_hist : input float

Other Parameters
----------------
nmp : input int, optional
Default: shape(xyz,1)
np : input int, optional
Default: len(id0_p)
nm : input int, optional
Default: len(id0_m)

Returns
-------
dist2_m : rank-1 array('d') with bounds (np*nm)
idx_m : int
'''

