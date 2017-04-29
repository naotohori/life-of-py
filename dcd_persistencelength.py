#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Created on 2013/03/12
@author: Naoto Hori
'''

import sys
from cafysis.file_io.dcd import DcdFile
import numpy as np
from scipy.optimize import curve_fit

if not len(sys.argv) in (6,):
    print ('Usage: % SCRIPT [input DCD] [skip frame] [first particle ID for calc] [gap] [output]')
    sys.exit(2)
    
frame_skip = int(sys.argv[2])
offset = int(sys.argv[3])
gap = int(sys.argv[4])
filepath_out = sys.argv[5]

dcd = DcdFile(sys.argv[1])
dcd.open_to_read()
dcd.read_header()

nmp = dcd.get_header().nmp_real

unit_len_sq = 0.0
n_unit_len_sq = 0
sum_cos_theta = [0.0] * nmp
num_n = [0] * nmp

dcd.skip(frame_skip)

first = True
while dcd.has_more_data() :
    data = dcd.read_onestep()
        
    for i in xrange(offset, nmp-1, gap) :
        for j in xrange(i+gap, nmp, gap) :
            if j == (i+gap):
                unit_len_sq += ( (data[offset-1][0] - data[offset-1+gap][0]) ** 2
                            +   +(data[offset-1][1] - data[offset-1+gap][1]) ** 2
                            +   +(data[offset-1][2] - data[offset-1+gap][2]) ** 2 )
                n_unit_len_sq += 1

            # i,j
            ij = (  data[i][0] * data[j][0]
                  + data[i][1] * data[j][1]
                  + data[i][2] * data[j][2] ) 
            # i-1,j
            i_1_j = (  data[i-1][0] * data[j][0]
                     + data[i-1][1] * data[j][1]
                     + data[i-1][2] * data[j][2] ) 
            # i,j-1
            i_j_1 = (  data[i][0] * data[j-1][0]
                     + data[i][1] * data[j-1][1]
                     + data[i][2] * data[j-1][2] ) 
            # i-1,j-1
            i_1_j_1 = (  data[i-1][0] * data[j-1][0]
                       + data[i-1][1] * data[j-1][1]
                       + data[i-1][2] * data[j-1][2] ) 

            cos_theta = ( ij - i_1_j - i_j_1 + i_1_j_1 ) 
            n = (j-i) / gap 
            sum_cos_theta[n] += cos_theta
            num_n[n] += 1
            
## Calculate average unit length
unit_len_sq = unit_len_sq / float(n_unit_len_sq)
unit_len = np.sqrt(unit_len_sq)

## Calculate average correlation
ij = []
cor = []
for i in xrange(nmp/gap):
    ij.append(float(i))
    if num_n[i] == 0:
        cor.append(1.0)
    else:
        cor.append( sum_cos_theta[i] / (float(num_n[i])*unit_len_sq) )

## Fitting
def func_exp(x,Lp):
    return np.exp(-(x*unit_len)/Lp)

para, dev = curve_fit(func_exp, ij, cor)

print ('Persistent length: %f nm' % (para[0]*0.1,))

## Output
f_out = open(filepath_out, 'w')

f_out.write('# L (unit length): %f\n' % unit_len)
f_out.write('# Lp (persistent length): %f\n' % para[0])
f_out.write('# pcov (fitting quality): %f\n' % dev[0])
f_out.write('#\n')
f_out.write('#  n   <cos>   exp(-n*L/Lp)\n')

for i in range(nmp/gap):
    f_out.write('%f %f %f %i\n' % (i, cor[i], np.exp(-(i*unit_len)/para[0]), num_n[i]) )

f_out.close()
