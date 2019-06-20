#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2011/03/17
@author: Naoto Hori
'''

import sys
from cafysis.file_io.dcd import DcdFile
from numpy import zeros, int32

if not len(sys.argv) in (4,5):
    print('Usage: % SCRIPT [input DCD] [cutoff] [output]')
    print('       % SCRIPT [input DCD] [cutoff] [output] [output PNG file (optional)]')
    sys.exit(2)
    
cutoff = float(sys.argv[2])
filepath_out = sys.argv[3]

dcd = DcdFile(sys.argv[1])
dcd.open_to_read()
dcd.read_header()

nmp = dcd.get_header().nmp_real
hist = zeros((nmp,nmp), dtype=int32)
cut_square = cutoff * cutoff

while dcd.has_more_data() :
    data = dcd.read_onestep()
    
    for i in range(nmp) :
        for j in range(i+1,nmp) :
            distance = ( (data[i][0] - data[j][0]) ** 2
                        +(data[i][1] - data[j][1]) ** 2
                        +(data[i][2] - data[j][2]) ** 2 )
            if distance <= cut_square :
                hist[i][j] += 1

file_out = open(filepath_out,'w')
for i in range(nmp) :
    for j in range(i+1,nmp) :
        file_out.write('%8i %8i %20i\n'%(i+1,j+1,hist[i][j]))
    file_out.write('\n')
file_out.close()    

if len(sys.argv) == 5:
    filepath_png = sys.argv[4]
    
    
    
    
        