#!/usr/bin/env python
#vim:fileencoding=UTF-8
'''
Created on 2016/08/05
@author: Naoto Hori

Calculate the centroid for particular part of molecule (defined by ID1 and ID2)
translate all coordinates so that the centroid comes to the origin (0,0,0).

Particles outside the box will be wrapped into the box in the origin.
'''

from cafysis.mtx_coord_transform import mtx_crd_transform
from cafysis.file_io.dcd import DcdFile
import sys
import math
import copy

if len(sys.argv) != 6:
    print('Usage: SCRIPT [input DCD] [ID domain begin] [ID domain end] [Box size] [output DCD]')
    sys.exit(2)

ID_DOM_INI = int(sys.argv[2]) - 1  # 重心を求める際に必要
ID_DOM_END = int(sys.argv[3]) - 1
BOXSIZE = float(sys.argv[4])

BOXMAX = 0.5 * BOXSIZE
BOXMIN = -0.5 * BOXSIZE

dcd = DcdFile(sys.argv[1])
dcd.open_to_read()
dcd_out = DcdFile(sys.argv[-1])
dcd_out.open_to_write()

# header
dcd.read_header()
header = dcd.get_header()
dcd_out.set_header(header)
dcd_out.write_header()

nmp = header.nmp_real

def calc_com_PBC(d):
    '''
    See https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions
    '''
    cos = [0.0] * 3
    sin = [0.0] * 3
    for xyz in d:
        for i in range(3):
            theta = (xyz[i] + BOXMAX) / BOXSIZE * 2 * math.pi
            cos[i] += math.cos(theta)
            sin[i] += math.sin(theta)

    com = [0.0] * 3
    for i in range(3):
        cos[i] = cos[i] / float(len(d))
        sin[i] = sin[i] / float(len(d))
        theta = math.atan2(-sin[i],-cos[i]) + math.pi
        com[i] = 0.5 * BOXSIZE * theta / math.pi - BOXMAX

    return com

def wrap(d):
    for i in range(len(d)):
        for j in range(3):
            p = d[i][j]
            if p > BOXMAX:
                d[i][j] = p - BOXSIZE * (int((p - BOXMAX)/BOXSIZE) + 1)
            elif p < BOXMIN:
                d[i][j] = p + BOXSIZE * (int((BOXMIN - p)/BOXSIZE) + 1)

while dcd.has_more_data() :
    data = dcd.read_onestep()
    
    ##########################################################
    #重心が原点に重なるように並進
    com = calc_com_PBC(data[ID_DOM_INI:ID_DOM_END+1])
    mtx_move = mtx_crd_transform()
    #mtx_move.reset()
    mtx_move.translation(-com[0],-com[1],-com[2])
   
    #for i in xrange(nmp):
    #    data[i][0:3] = trans_com.do_to_array(data[i])
    mtx_move.do_to_data(data)

    ##########################################################
    # Wrap particles outside the box at the orgin
    wrap(data)

    dcd_out.write_onestep(data) 
    
dcd.close()
dcd_out.close()
