#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@author: Naoto Hori
'''

import sys, os
from cafysis.file_io.dcd import DcdHeader,DcdFile
from cafysis.file_io.coord import CoordFile
from cafysis.file_io.psf import PsfFile
from cafysis.elements.psf import Atom, Psf

if not len(sys.argv) in (4,5):
    print 'Usage: % SCRIPT [input Coord] [#atom] [output DCD]'
    print '  or : % SCRIPT [input Coord] [#atom] [output DCD] [output PSF]'
    sys.exit(2)

filepath_crd = sys.argv[1]
nmp = int(sys.argv[2])
filepath_dcd = sys.argv[3]
filepath_psf = sys.argv[-1]

# PSF file
if len(sys.argv) == 5:
    # Generate PSF data
    psfdata = Psf()
    ires = 1
    for i in xrange(1, nmp+1):
        a = Atom()
        a.atom_id = i
        a.seg_name = 'R'
        if i%3 == 1:
            a.atom_name = 'C'
            a.atom_type = 'S'
        elif i%3 == 2:
            a.atom_name = 'N'
            a.atom_type = 'B'
        elif i%3 == 0:
            ires = ires + 1
            a.atom_name = 'O'
            a.atom_type = 'P'
            a.charge = -1.0
        a.res_id = ires
        a.res_name = 'R'
        psfdata.atoms.append(a)

    for i in xrange(2, nmp+1):
        if i%3 == 1: #S
            psfdata.bonds.append((i-1,i))
        elif i%3 == 2: #B
            psfdata.bonds.append((i-1,i))
        elif i%3 == 0: #P
            psfdata.bonds.append((i-2,i))

    # Write to file
    fpsf = PsfFile(filepath_psf)
    fpsf.open_to_write()
    fpsf.write_all(psfdata)
    fpsf.close()

# Calculate number of steps from file sieze
crdsize = os.path.getsize(filepath_crd)
num_step = crdsize / 8 / 3 / (nmp+1)

crd = CoordFile(filepath_crd, nmp)
crd.open_to_read()

header = DcdHeader()
header.nset = num_step
header.istart = 1
header.nstep_save = 1
header.nstep = 0
header.nunit_real = 0
header.delta = 0.05
header.title = ('','')
header.tempk = 0.0
header.lunit2mp = []
header.nmp_real = nmp

dcd = DcdFile(filepath_dcd)
dcd.open_to_write()
dcd.set_header(header)
dcd.write_header()

while crd.has_more_data():
    dcd.write_onestep(crd.read_onestep())

dcd.close()
crd.close()

