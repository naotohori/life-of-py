#!/usr/bin/env python
#vim:fileencoding=UTF-8
'''
Created on 2017/10/13
@author: Naoto Hori

Translate all coordinates so that the entire molecule is inside the box.
To do so, firstly find most lateral coordinates.

Particles outside the box will be wrapped into the box in the origin.
'''

from lop.mtx_coord_transform import mtx_crd_transform
from lop.file_io.dcd import DcdFile
import sys
import math
import copy

## Threshold of distance of neighboring beads
MAXD = 50.0

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

def find_max_min_PBC(d):

    max_xyz = d[0][:]
    min_xyz = d[0][:]
    pre = d[0][:]
    add = [0.0] * 3
    for xyz in d:
        for i in range(3):
            x = xyz[i] + add[i]
            if x - pre[i] > MAXD:
                x      += - BOXSIZE
                add[i] += - BOXSIZE
            elif x - pre[i] < -MAXD:
                x      +=   BOXSIZE
                add[i] +=   BOXSIZE
            pre[i] = x

            if max_xyz[i] < x:
                max_xyz[i] = x
            if x < min_xyz[i]:
                min_xyz[i] = x

    L = [0.0] * 3
    for i in range(3):
        L[i] = max_xyz[i] - min_xyz[i]

    return max_xyz, min_xyz, L

def wrap(d):
    for i in range(len(d)):
        for j in range(3):
            p = d[i][j]
            if p > BOXMAX:
                d[i][j] = p - BOXSIZE * (int((p - BOXMAX)/BOXSIZE) + 1)
            elif p < BOXMIN:
                d[i][j] = p + BOXSIZE * (int((BOXMIN - p)/BOXSIZE) + 1)

iframe = 0
while dcd.has_more_data() :
    data = dcd.read_onestep()
    
    ##########################################################
    # Find most lateral coordinates
    max_xyz, min_xyz, L = find_max_min_PBC(data[ID_DOM_INI:ID_DOM_END+1])

    if L[0] > BOXSIZE or L[1] > BOXSIZE or L[2] > BOXSIZE:
        print(('Warning: L exceeds BOXSIZE at frame %i' % iframe))

    mtx_move = mtx_crd_transform()
    #mtx_move.reset()
    mtx_move.translation(0.5*L[0] - max_xyz[0],
                         0.5*L[1] - max_xyz[1],
                         0.5*L[2] - max_xyz[2])
    
    #for i in xrange(nmp):
    #    data[i][0:3] = trans_com.do_to_array(data[i])
    mtx_move.do_to_data(data)

    ##########################################################
    # Wrap particles outside the box at the orgin
    wrap(data)

    # Copy the unit cell information
    dcd_out._header.unit_cell_xyz = dcd._header.unit_cell_xyz
    dcd_out._header.unit_cell_abc = dcd._header.unit_cell_abc

    dcd_out.write_onestep(data) 

    iframe += 1
    
dcd.close()
dcd_out.close()
