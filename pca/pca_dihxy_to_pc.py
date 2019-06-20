#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@author: Naoto Hori
'''

import sys
import struct
import os
import numpy as np

if len(sys.argv) < 5:
    print('Usage: % SCRIPT [input dcd] [ev file] [,ev file ....] [output average file] [output PC file]')
    sys.exit(2)
    
f_dihxy = open(sys.argv[1], 'rb')

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
    ev.append(np.array(ev_tmp))
    f_ev.close()
    
# Check for ev
num_dimension = len(ev[0])
for v in ev :
    if len(v) != num_dimension :
        print('len(v) != num_dimension, %i' % num_dimension)
        sys.exit(2)

# Calculate average structure (data_ave)
data_ave = [0.0 for i in range(num_dimension)]
num_data = 0

flg_more_data = True
while (flg_more_data):

    num = struct.unpack('i', f_dihxy.read(4))[0]

    for i in range(num):
        #x, y = struct.unpack('dd', f_dihxy.read(16))
        x, y = struct.unpack('ff', f_dihxy.read(8))
        data_ave[2*i]   += x
        data_ave[2*i+1] += y

    num_data += 1

    char = f_dihxy.read(4)
    if not char:
        flg_more_data = False
    else:
        f_dihxy.seek(-4,os.SEEK_CUR)
        flg_more_data = True

for i in range(len(data_ave)) :
    data_ave[i] /= float(num_data)
    f_ave.write('%12.5f\n' % data_ave[i])

f_ave.close()
            
f_dihxy.close()

f_dihxy = open(sys.argv[1], 'rb')

# Calculate principal coordinate
flg_more_data = True
while (flg_more_data):

    num = struct.unpack('i', f_dihxy.read(4))[0]

    data = []
    for i in range(num):
        #x, y = struct.unpack('dd', f_dihxy.read(16))
        x, y = struct.unpack('ff', f_dihxy.read(8))
        data.append(x - data_ave[2*i])
        data.append(y - data_ave[2*i+1])

    for v in ev:
        x = np.dot(v, np.array(data))
        f_out.write('%12.5f ' % x)
    f_out.write('\n')

    char = f_dihxy.read(4)
    if not char:
        flg_more_data = False
    else:
        f_dihxy.seek(-4,os.SEEK_CUR)
        flg_more_data = True

f_dihxy.close()
f_out.close()

