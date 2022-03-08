#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2017/03/21
@author: Naoto Hori
'''

import sys
import struct
from lop.torsion import torsion_xy
from lop.file_io.psf import PsfFile
from lop.file_io.dcd import DcdFile

if len(sys.argv) != 4:
    print ('\n Usage: SCRIPT [input DCD] [input PSF] [output]\n')
    sys.exit(2)

psffile = PsfFile(sys.argv[2])
psffile.open_to_read()
psf = psffile.read_all()

dcd = DcdFile(sys.argv[1])

f_out = open(sys.argv[-1],'wb')

id0_bb  = []
for i,a in enumerate(psf.atoms):
    if a.atom_name.strip() == 'P' or a.atom_name.strip() == 'S':
        id0_bb.append(i)

nbb = len(id0_bb)
ndih = nbb - 3

dih_mps = []
for idih in range(ndih):
    dih_mps.append( (id0_bb[idih], id0_bb[idih+1], id0_bb[idih+2], id0_bb[idih+3]) )

dcd.open_to_read()
dcd.read_header()

nstep = 0
while dcd.has_more_data() :
    nstep += 1

    data = dcd.read_onestep_np()

    f_out.write( struct.pack('i',ndih) )  # to make sure the data boundaries
    for mps in dih_mps:
        xy = torsion_xy(data[mps[0]], data[mps[1]], data[mps[2]], data[mps[3]])
        #f_out.write( struct.pack('dd',xy[0],xy[1]) )
        f_out.write( struct.pack('ff',xy[0],xy[1]) )

dcd.close()
