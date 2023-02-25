#!/usr/bin/env python
'''
Created on 2018/11/19
based on dcd_traj_distance.py
@author: Naoto Hori
'''

import sys
import math
from lop.file_io.dcd import DcdFile

if len(sys.argv) != 8:
    print(' Usage: % SCRIPT [input DCD] [boxsize x y z] [ID1] [ID2] [output] ')
    sys.exit(2)
    
dcd = DcdFile(sys.argv[1])
box = [float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])]
id1 = int(sys.argv[5]) - 1
id2 = int(sys.argv[6]) - 1

dcd.open_to_read()
dcd.read_header()

half = [0.5 * x for x in box]

def pb_neighbor(v):
    for i in range(3):
        if v[i] > half[i]:
            v[i] = v[i] - box[i]
        elif v[i] < -half[i]:
            v[i] = v[i] + box[i]
    return v

f_out = open(sys.argv[-1],'w')
#nframe = 0
while dcd.has_more_data() :
    data = dcd.read_onestep()
    #nframe += 1
    v = [data[id1][0] - data[id2][0], 
         data[id1][1] - data[id2][1],
         data[id1][2] - data[id2][2] ]

    print(v)
    v = pb_neighbor(v)

    d = math.sqrt( v[0]**2 + v[1]**2 + v[2]**2 )
    print(v, d)

    f_out.write('%.2f\n' % d)
dcd.close()

f_out.close()
