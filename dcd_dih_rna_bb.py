#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2017/03/21
@author: Naoto Hori
'''

import sys
from cafysis.torsion import torsion
from cafysis.file_io.psf import PsfFile
from cafysis.file_io.dcd import DcdFile

if len(sys.argv) != 4:
    print ('\n Usage: SCRIPT [input DCD] [input PSF] [output]\n')
    sys.exit(2)

psffile = PsfFile(sys.argv[2])
psffile.open_to_read()
psf = psffile.read_all()

dcd = DcdFile(sys.argv[1])

f_out = open(sys.argv[-1],'w')

id0_bb  = []
for i,a in enumerate(psf.atoms):
    if a.atom_name.strip() == 'P' or a.atom_name.strip() == 'S':
        id0_bb.append(i)

nbb = len(id0_bb)

dih_mps = []
for idih in range(nbb - 3):
    dih_mps.append( (id0_bb[idih], id0_bb[idih+1], id0_bb[idih+2], id0_bb[idih+3]) )

dcd.open_to_read()
dcd.read_header()

nstep = 0
while dcd.has_more_data() :
    nstep += 1

    data = dcd.read_onestep_np()

    for mps in dih_mps:
        f_out.write(' %8.3f' %  torsion(data[mps[0]], data[mps[1]], data[mps[2]], data[mps[3]], flg_degree=True, flg_360=False) )

    f_out.write('\n')

dcd.close()
