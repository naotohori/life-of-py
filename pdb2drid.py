#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Created on 2015/05/26
@author: Naoto Hori
'''

import sys
from cafysis.file_io.pdb import PdbFile
from cafysis.file_io.drid import DridFile, DridHeader
from cafysis.lib_f2py import py_drid
import numpy as np

if len(sys.argv) != 5:
    print 'Usage: % SCRIPT [input PDB] [centroid set file] [atom set file] [output DRID]'
    sys.exit(2)

pdb = PdbFile(sys.argv[1])
pdb.open_to_read()
chains = pdb.read_all()
pdb.close()

centroids = []
for l in open(sys.argv[2],'r'):
    centroids.append(int(l))

atoms = []
for l in open(sys.argv[3],'r'):
    atoms.append(int(l))

bonds = []

drid = DridFile(sys.argv[-1])
drid.open_to_write()

h = DridHeader()
h.title = 'ver1'
h.centroids = centroids
h.atoms = atoms
h.bonds = bonds
drid.set_header(h)
drid.write_header()

nmp = 0
for c in chains:
    nmp += c.num_atom()

data = np.zeros((nmp,3))
imp = 0
for c in chains:
    for i in range(c.num_atom()):
        a = c.get_atom(i)
        data[imp,0] = a.xyz.x
        data[imp,1] = a.xyz.y
        data[imp,2] = a.xyz.z
        imp += 1

'''
F2PY
munuxi = drid(x,centroids,atoms,nx=shape(x,1),nc=len(centroids),na=len(atoms))
'''
munuxi = py_drid.drid( data.T, centroids, atoms) 

drid.write_onestep( munuxi )

drid.close()
