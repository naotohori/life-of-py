#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2016/10/13
Based on dcd_traj_SaltDistribution_nt.py

@author: Naoto Hori
'''

import sys
import math
import numpy as np
from cafysis.lib_f2py import py_distance_count_within
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

r_kT = 7.2684

n_sum_Mg = np.zeros(nP)
n_sum_K  = np.zeros(nP)
n_sum_Cl = np.zeros(nP)
nsq_sum_Mg = np.zeros(nP)
nsq_sum_K  = np.zeros(nP)
nsq_sum_Cl = np.zeros(nP)

dcd.open_to_read()
dcd.read_header()

nstep = 0
while dcd.has_more_data() :
    nstep += 1

    if nstep % 10000 == 0:
        print nstep

    data = dcd.read_onestep_npF()

    if len(id0_Mg) > 0:
        num = py_distance_count_within.distance_count_within(data, id0_P, id0_Mg, r_kT)
        n_sum_Mg += num
        nsq_sum_Mg += num ** 2

    if len(id0_K) > 0:
        num = py_distance_count_within.distance_count_within(data, id0_P, id0_K, r_kT)
        n_sum_K += num
        nsq_sum_K += num ** 2

    if len(id0_Cl) > 0:
        num = py_distance_count_within.distance_count_within(data, id0_P, id0_Cl, r_kT)
        n_sum_Cl += num
        nsq_sum_Cl += num ** 2

dcd.close()

f_out.write('#nstep: %i\n#\n' % (nstep,))

div_factor = float(nstep) * (4.0 / 3.0 * math.pi * (r_kT**3) * 6.02214 * 0.0001)
div2_factor = float(nstep) * (4.0 / 3.0 * math.pi * (r_kT**3) * 6.02214 * 0.0001) **2

for iP in range(nP):

    ave_K  = n_sum_K[iP] / div_factor 
    dev_K  = math.sqrt( nsq_sum_K[iP] / div2_factor - ave_K ** 2)
    ave_Cl = n_sum_Cl[iP] / div_factor
    dev_Cl = math.sqrt( nsq_sum_Cl[iP] / div2_factor - ave_Cl ** 2)

    if len(id0_Mg) > 0 and len(id0_K) > 0:
        ave_Mg = n_sum_Mg[iP] / div_factor
        dev_Mg = math.sqrt( nsq_sum_Mg[iP] / div2_factor - ave_Mg ** 2)

        f_out.write('%5i %g %g %g %g %g %g\n' % (iP+1,ave_Mg,ave_K,ave_Cl,dev_Mg,dev_K,dev_Cl))

    elif len(id0_K) > 0:
        f_out.write('%5i %g %g %g %g %g %g\n' % (iP+1,0.0,ave_K,ave_Cl,0.0,dev_K,dev_Cl))

f_out.close()

