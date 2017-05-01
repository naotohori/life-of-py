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
sum_cos_theta = [0.0] * (nmp-1)
num_n = [0] * (nmp-1)

dcd.skip(frame_skip)

first = True
while dcd.has_more_data() :
    data = dcd.read_onestep_np()
        
    for i in xrange(offset, nmp-1, gap) :

        vi = data[i] - data[i-1]
        unit_len_sq += np.dot(vi,vi)
        n_unit_len_sq += 1

        for j in xrange(i+gap, nmp, gap) :

            vj = data[j] - data[j-1]

            cos_theta = np.dot(vi, vj)

            n = (j-i) / gap 
            sum_cos_theta[n] += cos_theta
            num_n[n] += 1
            
## Calculate average unit length
unit_len_sq = unit_len_sq / float(n_unit_len_sq)
unit_len = np.sqrt(unit_len_sq)

## Calculate average correlation
ij = []
cor = []
for i in xrange(nmp/gap-1):
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

for i in range(nmp/gap-1):
    f_out.write('%f %f %f %i\n' % (i, cor[i], np.exp(-(i*unit_len)/para[0]), num_n[i]) )

f_out.close()
