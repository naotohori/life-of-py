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

if len(sys.argv) not in (4,5):
    print('Usage: % SCRIPT [input DCD] [mask file] [output DRID]')
    print('  or : % SCRIPT [input DCD] [mask file] [# solute in DCD] [output DRID]')
    sys.exit(2)

dcd = DcdFile(sys.argv[1])
dcd.open_to_read()
dcd.read_header()

mask = []
for l in open(sys.argv[2]):
    lsp = l.strip().split()
    a = []
    for x in lsp:
        a.append( float(x) ) 
    mask.append( a )

if len(sys.argv) == 5:
    flg_solute = True
    nsolute = int(sys.argv[3])
    if len(mask) != nsolute or len(mask[0]) != nsolute:
        print('Error: # solute and size of mask is not consistent')
        sys.exit(2)

drid = DridFile(sys.argv[-1])
drid.open_to_write()

h = DridHeader()
h.title = 'ver2'
h.version = 2
h.mask = np.array( mask )
drid.set_header(h)
drid.write_header()

#iframe = 0
while dcd.has_more_data() :

    if flg_solute:
        data = dcd.read_onestep_np_solute(nsolute)
    else:
        data = dcd.read_onestep_np()
    
    '''
    F2PY
    '''
    munuxi = py_drid.drid( data.T, h.mask.T )

    drid.write_onestep( munuxi )

    '''
    PYTHON version (time consuming; for debugging)
    '''
    '''
    nc = h.mask.shape[0]
    mu = [0.0]*nc
    nsum = [0]*nc
    for ic in range(nc):
        for ia in range(nc):
            if h.mask[ic, ia] < 0:
                pass
            else:
                d = math.sqrt( np.square( (data[ia,:] - data[ic,:])).sum() )
                mu[ic] += 1.0 / d
                nsum[ic] += 1

    for ic in range(nc):
        mu[ic] = mu[ic] / float(nsum[ic])

    nu = [0.0]*nc
    xi = [0.0]*nc
    for ic in range(nc):
        for ia in range(nc):
            if h.mask[ic,ia] < 0:
                pass
            else:
                d = math.sqrt( np.square( (data[ia,:] - data[ic,:])).sum() )
                dd = 1.0/d - mu[ic]
                nu[ic] += dd ** 2
                xi[ic] += dd ** 3

    for ic in range(nc):
        nu[ic] = math.sqrt( nu[ic] / float(nsum[ic]) )
        if xi[ic] < 0.0:
            xi[ic] = - ( (abs(xi[ic]) / float(nsum[ic])) ** (1./3.) )
        else:
            xi[ic] = (xi[ic] / float(nsum[ic])) ** (1./3.) 
        #xi[ic] = (xi[ic] / float(nsum[ic]) + 0j) ** (1./3.)     
        ##!!!!! x ** (1./3.) causes unstability when x < 0

    # Check Fortran vs Python
    for ic in range(nc):
        print munuxi[ic], mu[ic], munuxi[ic]-mu[ic]
    for ic in range(nc):
        print munuxi[nc+ic], nu[ic], munuxi[nc+ic]-nu[ic]
    for ic in range(nc):
        print munuxi[2*nc+ic], xi[ic], munuxi[2*nc+ic]-xi[ic]
    '''

    #drid.write_onestep( np.concatenate((mu,nu,xi)) )

    #print iframe
    #iframe += 1
