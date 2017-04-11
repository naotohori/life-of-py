#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Created on 2015/05/15
@author: Naoto Hori
'''

import sys
import math
from cafysis.file_io.dcd import DcdFile
from cafysis.file_io.drid import DridFile, DridHeader
import cafysis.lib_f2py.py_drid as py_drid
import numpy as np

if len(sys.argv) != 5:
    print 'Usage: % SCRIPT [input DCD] [centroid set file] [atom set file] [output DRID]'
    sys.exit(2)

dcd = DcdFile(sys.argv[1])
dcd.open_to_read()
dcd.read_header()

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

#nc = len(centroids)
#na = len(atoms)
#mu = []
#iframe = 0
while dcd.has_more_data() :
    data = dcd.read_onestep_np()
    
    '''
    F2PY
    munuxi = drid(x,centroids,atoms,nx=shape(x,1),nc=len(centroids),na=len(atoms))
    '''
    munuxi = py_drid.drid( data.T, centroids, atoms) 
    drid.write_onestep( munuxi )

    '''
    # PYTHON version (time consuming)
    d_inv = np.ma.empty((nc, na)) # [1:nc, 1:na]
    d_inv.mask = False

    nsum = []
    for ic,c in enumerate(centroids):
        coord_c = data[c-1,:]

        ## add here treatment for bonds
        exclude = [c,]

        n = 0
        for ia, a in enumerate(atoms):
            if a in exclude:
                d_inv[ic,ia] = np.ma.masked
            else:
                d_inv[ic,ia] = float(1) / math.sqrt( np.square((data[a-1,:] - coord_c)).sum() )
                n += 1
        nsum.append(float(n))

    mu = d_inv.sum(axis=1)  # [1:nc]
    mu = mu / nsum

    dd = np.ma.empty((nc,na))
    dd.mask = d_inv.mask
    for ic in range(nc):
        dd[ic,:] = d_inv[ic,:] - mu[ic]

    nu = np.sqrt( np.square(dd).sum(axis=1) / nsum )
    xi = np.power( np.power(dd,3).sum(axis=1) / nsum, float(1)/3)

    drid.write_onestep( np.concatenate((mu,nu,xi)) )
    '''
    #print iframe
    #iframe += 1
