#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Created on 2013/03/12
@author: Naoto Hori
'''

import sys
import argparse
from cafysis.file_io.dcd import DcdFile
import numpy as np
from scipy.optimize import curve_fit

parser = argparse.ArgumentParser(description='Calculate persistence length from DCD trajectory',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--gap', dest='gap', default=1,
                    action='store', type=int, help='Gap')
parser.add_argument('--skip', dest='frame_skip', default=0,
                    action='store', type=int, help='number of frames to be skipped')
parser.add_argument('--first', dest='offset', default=0,
                    action='store', type=int, help='first particle ID')
parser.add_argument('--PosCor', dest='flg_poscor', default=False,
                    action='store_true', help='Consider only the length which has positive correlation')

parser.add_argument('dcd', help='target DCD file')
parser.add_argument('out', help='output file')

args = parser.parse_args()

################################################################################

dcd = DcdFile(args.dcd)
dcd.open_to_read()
dcd.read_header()

nmp = dcd.get_header().nmp_real

unit_len_sq = 0.0
n_unit_len_sq = 0
sum_cos_theta = [0.0] * (nmp-1)
num_n = [0] * (nmp-1)

dcd.skip(args.frame_skip)

first = True
while dcd.has_more_data() :
    data = dcd.read_onestep_np()
        
    #print 'Start i loop'
    for i in xrange(args.offset+args.gap, nmp-args.gap, args.gap) :

        #print i, i-args.gap
        vi = data[i] - data[i-args.gap]
        unit_len_sq += np.dot(vi,vi)
        n_unit_len_sq += 1

        #print 'Start j loop'
        for j in xrange(i+args.gap, nmp, args.gap) :

            #print '    ',j, j-args.gap
            vj = data[j] - data[j-args.gap]

            cos_theta = np.dot(vi, vj)

            n = (j-i) / args.gap 
            sum_cos_theta[n] += cos_theta
            num_n[n] += 1
            
## Calculate average unit length
unit_len_sq = unit_len_sq / float(n_unit_len_sq)
unit_len = np.sqrt(unit_len_sq)

## Calculate average correlation
cor = []
for i in xrange(nmp/args.gap-1):
    if num_n[i] == 0:
        cor.append(1.0)
    else:
        cor.append( sum_cos_theta[i] / (float(num_n[i])*unit_len_sq) )

if args.flg_poscor:
    n_cor = 0
    ij = []
    for i, c in enumerate(cor):
        if c < 0.0:
            break
        n_cor += 1
        ij.append(float(i))
else:
    n_cor = len(cor)
    ij = [float(i) for i in range(n_cor)]

## Fitting
def func_exp(x,Lp):
    return np.exp(-(x*unit_len)/Lp)

para, dev = curve_fit(func_exp, ij, cor[:n_cor])

print ('Persistent length: %f nm' % (para[0]*0.1,))

## Output
f_out = open(args.out, 'w')

f_out.write('# L (unit length): %f\n' % unit_len)
f_out.write('# Lp (persistence length): %f\n' % para[0])
f_out.write('# pcov (fitting quality): %f\n' % dev[0])
f_out.write('#\n')
f_out.write('#  n   <cos>   exp(-n*L/Lp)\n')

for i in range(nmp/args.gap-1):
    f_out.write('%5i %6.2f %7.4f %7.4f %12i\n' % (i, i*unit_len, cor[i], np.exp(-(i*unit_len)/para[0]), num_n[i]) )

f_out.close()
