#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2016/10/04
@author: Naoto Hori

Ref. pdb_centroid_origin_PBC.py
'''

import sys
import math
import numpy as np
from cafysis.file_io.dcd import DcdFile
from cafysis.file_io.psf import PsfFile
from cafysis.mtx_coord_transform import mtx_crd_transform

if len(sys.argv) != 5:
    print ('\n Usage: SCRIPT [input DCD] [input PSF] [Box size] [output]\n')
    sys.exit(2)
    
dcd = DcdFile(sys.argv[1])

psffile = PsfFile(sys.argv[2])
psffile.open_to_read()
psf = psffile.read_all()

BOXSIZE = float(sys.argv[3])
BOXMAX = 0.5 * BOXSIZE
BOXMIN = -0.5 * BOXSIZE

def calc_com_PBC(d):
    '''
    See https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions
    '''
    cos = [0.0] * 3
    sin = [0.0] * 3
    for xyz in d:
        for i in range(3):
            theta = (xyz[i] + BOXMAX) / BOXSIZE * 2 * math.pi
            cos[i] += math.cos(theta)
            sin[i] += math.sin(theta)

    com = [0.0] * 3
    for i in range(3):
        cos[i] = cos[i] / float(len(d))
        sin[i] = sin[i] / float(len(d))
        theta = math.atan2(-sin[i],-cos[i]) + math.pi
        com[i] = 0.5 * BOXSIZE * theta / math.pi - BOXMAX

    return com

def wrap(d):
    for i in range(len(d)):
        for j in range(3):
            p = d[i][j]
            if p > BOXMAX:
                d[i][j] = p - BOXSIZE * (int((p - BOXMAX)/BOXSIZE) + 1)
            elif p < BOXMIN:
                d[i][j] = p + BOXSIZE * (int((BOXMIN - p)/BOXSIZE) + 1)

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

dr = 1.0
bins = [(x*dr)**2 for x in range(175+1)]
nbin = len(bins) - 1

hist_Mg = np.zeros((nbin,))
hist_K  = np.zeros((nbin,))
hist_Cl = np.zeros((nbin,))
hist_P  = np.zeros((nbin,))

dcd.open_to_read()
dcd.read_header()

nstep = 0
while dcd.has_more_data() :
    nstep += 1
    #if nstep % 1000 == 0:
        #print nstep

    data = dcd.read_onestep()

    #重心が原点に重なるように並進
    com = calc_com_PBC(data[0:587])
    mtx_move = mtx_crd_transform()
    #mtx_move.reset()
    mtx_move.translation(-com[0],-com[1],-com[2])
    mtx_move.do_to_data(data)

    # Wrap particles outside the box at the orgin
    wrap(data)

    if len(id0_Mg) > 0:
        darray = []
        for iMg in id0_Mg:
            d2 = data[iMg][0] ** 2 + data[iMg][1] ** 2 + data[iMg][2] ** 2
            darray.append(d2)

        hist, B = np.histogram(darray, bins=bins, density=False)
        hist_Mg[:] += hist[:]

    if len(id0_K) > 0:
        darray = []
        for iK in id0_K:
            d2 = data[iK][0] ** 2 + data[iK][1] ** 2 + data[iK][2] ** 2
            darray.append(d2)

        hist, B = np.histogram(darray, bins=bins, density=False)
        hist_K[:] += hist[:]

    if len(id0_P) > 0:
        darray = []
        for iP in id0_P:
            d2 = data[iP][0] ** 2 + data[iP][1] ** 2 + data[iP][2] ** 2
            darray.append(d2)

        hist, B = np.histogram(darray, bins=bins, density=False)
        hist_P[:] += hist[:]

    if len(id0_Cl) > 0:
        darray = []
        for iCl in id0_Cl:
            d2 = data[iCl][0] ** 2 + data[iCl][1] ** 2 + data[iCl][2] ** 2
            darray.append(d2)

        hist, B = np.histogram(darray, bins=bins, density=False)
        hist_Cl[:] += hist[:]

    #if nstep == 10000:
    #    break

dcd.close()

f_out.write('#nstep: %i\n' % (nstep,))
f_out.write('# 1, 2 : r1, r2\n')
f_out.write('# 3, 4, 5, 6 : n_P, n_Mg, n_K, n_Cl\n')
f_out.write('# 7, 8, 9,10 : c_P, c_Mg, c_K, c_Cl\n')
f_out.write('#11,12,13,14 : ex_P, ex_Mg, ex_K, ex_Cl\n')
f_out.write('#15,16,17,18 : sum_ex_P, sum_ex_Mg, sum_ex_K, sum_ex_Cl\n')

fn = float(nstep)
r1 = math.sqrt(bins[-6])
r2 = math.sqrt(bins[-1])
v_bulk = 4.0 / 3.0 * math.pi * (r2 ** 3 - r1 ** 3)
c_bulk_P  = 0
c_bulk_Mg = sum(hist_Mg[-5:]) / v_bulk / fn
c_bulk_K  = sum(hist_K[-5:])  / v_bulk / fn
c_bulk_Cl = sum(hist_Cl[-5:]) / v_bulk / fn
f_out.write('#\n')
f_out.write('#Bulk (r1,r2) = %f %f\n' % (r1,r2))
f_out.write('#c_bulk_Mg : %g\n' % c_bulk_Mg)
f_out.write('#c_bulk_K  : %g\n' % c_bulk_K)
f_out.write('#c_bulk_Cl : %g\n' % c_bulk_Cl)

sum_n_P = 0 ; sum_n_Mg = 0 ; sum_n_K = 0 ; sum_n_Cl = 0 
sum_ex_P = 0.0 ; sum_ex_Mg = 0.0 ; sum_ex_K = 0.0 ; sum_ex_Cl = 0.0 

for ib in range(nbin):
    n_P  = hist_P[ib]
    n_Mg = hist_Mg[ib]
    n_K  = hist_K[ib]
    n_Cl = hist_Cl[ib]
    sum_n_P  += n_P
    sum_n_Mg += n_Mg
    sum_n_K  += n_K
    sum_n_Cl += n_Cl

    ave_P  = n_P  / fn
    ave_Mg = n_Mg / fn
    ave_K  = n_K  / fn
    ave_Cl = n_Cl / fn

    r1 = math.sqrt(bins[ib])
    r2 = math.sqrt(bins[ib+1])
    v = 4.0 / 3.0 * math.pi * (r2 ** 3 - r1 ** 3)   # Angstrome

    c_P  = ave_P  / v
    c_Mg = ave_Mg / v
    c_K  = ave_K  / v
    c_Cl = ave_Cl / v

    ex_P  = c_P  - c_bulk_P
    ex_Mg = c_Mg - c_bulk_Mg
    ex_K  = c_K  - c_bulk_K
    ex_Cl = c_Cl - c_bulk_Cl

    #sum_ex_P  += ex_P  * v
    #sum_ex_Mg += ex_Mg * v
    #sum_ex_K  += ex_K  * v
    #sum_ex_Cl += ex_Cl * v
    v2 = 4.0 / 3.0 * math.pi * (r2 ** 3)   # Angstrome
    sum_ex_P  = sum_n_P  / fn - c_bulk_P  * v2
    sum_ex_Mg = sum_n_Mg / fn - c_bulk_Mg * v2
    sum_ex_K  = sum_n_K  / fn - c_bulk_K  * v2
    sum_ex_Cl = sum_n_Cl / fn - c_bulk_Cl * v2

    f_out.write('%4.0f %4.0f %i %i %i %i %g %g %g %g %g %g %g %g %g %g %g %g %i %i %i %i\n' %
                 (r1,r2,
                  n_P,n_Mg,n_K,n_Cl,
                  c_P,c_Mg,c_K,c_Cl,
                  ex_P,ex_Mg,ex_K,ex_Cl,
                  sum_ex_P,sum_ex_Mg,sum_ex_K,sum_ex_Cl,
                  sum_n_P, sum_n_Mg, sum_n_K, sum_n_Cl))
f_out.close()
