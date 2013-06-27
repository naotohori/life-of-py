#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2012/04/30
@author: Naoto Hori
'''

if __name__ == '__main__':
    pass

import sys
import math
import numpy as np

if len(sys.argv) != 9:
    print ("Usage: SCRIPT [ts file] [x min] [x max] [x bin] [y min] [y max] [y bin] [out]")
    sys.exit(2)
   
K_B = 1.98 * 0.001
TEMP = 300.0
KBT = K_B * TEMP
ene = []
pc1 = []
pc2 = []

pc1_min = float(sys.argv[2])
pc1_max = float(sys.argv[3])
pc1_wid = float(sys.argv[4])
pc2_min = float(sys.argv[5])
pc2_max = float(sys.argv[6])
pc2_wid = float(sys.argv[7])

num_out_range = 0
for line in open(sys.argv[1], 'r'):
    if line.find('#') != -1:
        continue
    linesp = line.split()
    e = float(linesp[3])
    x1 = float(linesp[7])
    x2 = float(linesp[8])
    if (x1 < pc1_min or x1 > pc1_max or x2 < pc2_min or x2 > pc2_max):
        num_out_range += 1
        print ("#Out of the range: %s" % line.strip())
        continue
    ene.append(float(linesp[3]))
    pc1.append(float(linesp[7]))
    pc2.append(float(linesp[8]))
    
print("### Number of out-of-range")
print ("%i" % num_out_range)
print("### Confirmation: number of data")
print("ene: %i" % len(ene))
print("pc1: %i" % len(pc1))
print("pc2: %i" % len(pc2))
print("### Confirmation: min and max")
print("ene: %f - %f" % (min(ene), max(ene)))
print("pc1: %f - %f" % (min(pc1), max(pc1)))
print("pc2: %f - %f" % (min(pc2), max(pc2)))

pc1_bins = np.arange(pc1_min, pc1_max, pc1_wid)
pc2_bins = np.arange(pc2_min, pc2_max, pc2_wid)
#print("##### pc1_bins")
#print(pc1_bins)
#print("##### pc2_bins")
#print(pc2_bins)

pc1_inds = np.digitize(pc1, pc1_bins)
pc2_inds = np.digitize(pc2, pc2_bins)
#print("##### pc1_inds")
#print(pc1_inds)
#print("##### pc2_inds")
#print(pc2_inds)

z = np.zeros((len(pc1_bins),len(pc2_bins)))
for i,e in enumerate(ene):
    z[pc1_inds[i]-1, pc2_inds[i]-1] += math.exp(-e/KBT)
    
f_out = open(sys.argv[-1],'w')
for i,x in enumerate(pc1_bins):
    for j,y in enumerate(pc2_bins):
        if z[i,j] <= 0.0:
            f = 0.0
        else:
            f = -KBT * math.log(z[i,j])
        f_out.write('%f %f %f\n' % (0.5*(2*x+pc1_wid), 0.5*(2*y+pc1_wid), f))
    f_out.write('\n')
f_out.close()
        
