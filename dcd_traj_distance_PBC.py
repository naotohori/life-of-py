#!/usr/bin/env python
'''
Created on 2018/11/19
based on dcd_traj_distance.py
@author: Naoto Hori
'''

import sys
import math
from cafysis.file_io.dcd import DcdFile

if len(sys.argv) != 6:
    print(' Usage: % SCRIPT [input DCD] [boxsize] [ID1] [ID2] [output PDB] ')
    sys.exit(2)
    
dcd = DcdFile(sys.argv[1])
boxsize = float(sys.argv[2])
id1 = int(sys.argv[3]) - 1
id2 = int(sys.argv[4]) - 1

dcd.open_to_read()
dcd.read_header()

box = [boxsize]*3
half = [0.5*boxsize]*3

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

    v = pb_neighbor(v)

    d = math.sqrt( v[0]**2 + v[1]**2 + v[2]**2 )

    f_out.write('%.2f\n' % d)
dcd.close()

f_out.close()
