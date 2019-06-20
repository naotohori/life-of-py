#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@author: Naoto Hori
'''

import sys
from cafysis.file_io.dcd import DcdFile

if len(sys.argv) < 5:
    print('Usage: % SCRIPT [input dcd] [ev file] [,ev file ....] [output average file] [output PC file]')
    sys.exit(2)
    
f_dcd = DcdFile(sys.argv[1])
f_dcd.open_to_read()
f_dcd.read_header()

f_out = open(sys.argv[-1], 'w')
f_ave = open(sys.argv[-2], 'w')

# Read eigen values
num_ev = len(sys.argv) - 4
ev = []
for i in range(num_ev) :
    f_ev = open(sys.argv[i+2], 'r')
    ev_tmp = [] 
    for line in f_ev :
        if line.find('#') != -1: continue
        ev_tmp.append(float(line.strip()))
    ev.append(ev_tmp)
    f_ev.close()
    
# Check for ev
num_dimension = len(ev[0])
for v in ev :
    if len(v) != num_dimension :
        print('len(v) != num_dimension, %i' % num_dimension)
        sys.exit(2)

#debug
#for i in xrange(num_ev):
#    print '# ev %i' % i
#    for value in ev[i] :
#        print value

# Calculate average structure (data_ave)
data_ave = [0.0 for i in range(num_dimension)]
num_data = 0
while f_dcd.has_more_data() :
    num_data += 1
    data = f_dcd.read_onestep()
    idx = 0
    for atom in data :
        for xyz in atom :
            data_ave[idx] += xyz
            idx += 1
for i in range(len(data_ave)) :
    data_ave[i] /= float(num_data)
    f_ave.write('%12.5f\n' % data_ave[i])

f_ave.close()
            
f_dcd.close()
f_dcd.open_to_read()
f_dcd.read_header()

# Calculate principal coordinate
while f_dcd.has_more_data() :
    data = f_dcd.read_onestep()
    
    for v in ev :
        x = 0.0
        idx = 0
        for atom in data :
            for xyz in atom :
                x += (xyz - data_ave[idx]) * v[idx]
                idx += 1
        f_out.write('%12.5f ' % x)
    f_out.write('\n')
    
f_dcd.close()
f_out.close()

