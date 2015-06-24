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
from cafysis.lib_f2py import py_drid
#import py_drid
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
    munuxi = drid(x,centroids,atoms,nx=shape(x,1),nc=len(centroids),na=len(atoms))
    '''
    munuxi = py_drid.drid( data.T, centroids, atoms) 

    drid.write_onestep( munuxi )
    #print iframe
    #iframe += 1
