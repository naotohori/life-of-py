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

if len(sys.argv) != 4:
    print 'Usage: % SCRIPT [input PDB] [mask file] [output DRID]'
    sys.exit(2)

pdb = PdbFile(sys.argv[1])
pdb.open_to_read()
chains = pdb.read_all()
pdb.close()

mask = []
for l in open(sys.argv[2]):
    lsp = l.strip().split()
    a = []
    for x in lsp:
        a.append( float(x) ) 
    mask.append( a )


drid = DridFile(sys.argv[-1])
drid.open_to_write()

h = DridHeader()
h.title = 'ver2'
h.version = 2
h.mask = np.array( mask )
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
#munuxi = py_drid.drid( data.T, centroids, atoms) 
munuxi = py_drid.drid( data.T, h.mask.T )

drid.write_onestep( munuxi )

drid.close()
