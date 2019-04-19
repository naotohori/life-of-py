#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2015/06/17
@author: Naoto Hori
'''
import sys
import math
import numpy as np
from cafysis.lib_f2py import py_dcd_r2_histogram
from cafysis.file_io.dcd import DcdFile

if len(sys.argv) != 5:
    print ('\n Usage: SCRIPT [input DCD] [nmp] [nskip] [max r]\n')
    sys.exit(2)

print('#',sys.argv)
    
nmp = int(sys.argv[2])
nskip = int(sys.argv[3])
max_r = int(sys.argv[4])
nbin = max_r * 10
bin_edges = [x*0.1 for x in range(nbin+1)]
bin_edges_sq = [(x*0.1)**2 for x in range(nbin+1)]

dcd = DcdFile(sys.argv[1])
dcd.open_to_read()
dcd.read_header()

hist_all = [0] * nbin
icount = 0
while dcd.has_more_data() :

    data = dcd.read_onestep_npF()

    hist = py_dcd_r2_histogram.dcd_r2_histogram(data[:,:nmp], bin_edges_sq)
   
    hist_all += hist
    icount += 1
    
    if dcd.has_more_data():
        dcd.skip(nskip-1)

dcd.close()

for i in range(nbin):
    print(bin_edges[i], bin_edges[i+1], hist_all[i] / float(icount))


'''
hist = dcd_r2_histogram(xyz,bin_edges,[nmp,nbin])

Wrapper for ``dcd_r2_histogram``.

Parameters
----------
xyz : input rank-2 array('d') with bounds (3,nmp)
bin_edges : input rank-1 array('d') with bounds (nbin + 1)

Other Parameters
----------------
nmp : input int, optional
    Default: shape(xyz,1)
nbin : input int, optional
    Default: (len(bin_edges)-1)

Returns
-------
hist : rank-1 array('i') with bounds (nbin)
'''
