#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2015/06/05
@author: Naoto Hori
'''
import sys
import math
import numpy as np
from cafysis.lib_f2py import py_count_bound
from cafysis.file_io.dcd import DcdFile
from cafysis.file_io.psf import PsfFile

if len(sys.argv) != 4:
    print ('\n Usage: SCRIPT [input DCD] [input PSF] [output prefix]\n')
    sys.exit(2)
    
psffile = PsfFile(sys.argv[2])
psffile.open_to_read()
psf = psffile.read_all()


id0_Mg = []
id0_P  = []
id0_K  = []
for i,a in enumerate(psf.atoms):
    if a.atom_name.strip() == 'P':
        id0_P.append(i)
    elif a.atom_name.strip() == 'Mg':
        id0_Mg.append(i)
    elif a.atom_name.strip() == 'K':
        id0_K.append(i)


f_out_Mg  = open(sys.argv[-1]+'.boundsalt.Mg','w')
f_out_K   = open(sys.argv[-1]+'.boundsalt.K','w')
f_out_total = open(sys.argv[-1]+'.boundsalt.total','w')

r_count_Mg = 0.5*(2.1+1.5852) + 1.5
r_count_K  = 0.5*(2.1+5.3160) + 1.5

dcd = DcdFile(sys.argv[1])
dcd.open_to_read()
dcd.read_header()

while dcd.has_more_data() :

    data = dcd.read_onestep_npF()

    if len(id0_Mg) > 0:
        count_Mg = py_count_bound.count_bound( data, id0_P, id0_Mg, r_count_Mg)

    for c in count_Mg:
        f_out_Mg.write('%i ' % (c,))
    f_out_Mg.write('\n')

    if len(id0_K) > 0:
        count_K = py_count_bound.count_bound( data, id0_P, id0_K, r_count_K)

    for c in count_K:
        f_out_K.write('%i ' % (c,))
    f_out_K.write('\n')

    f_out_total.write('%i %i\n' % (sum(count_Mg), sum(count_K),))

dcd.close()

f_out_Mg.close()
f_out_K.close()
f_out_total.close()
