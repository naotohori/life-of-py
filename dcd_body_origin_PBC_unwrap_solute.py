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

dcd = DcdFile(sys.argv[1])
dcd.open_to_read()
dcd_out = DcdFile(sys.argv[-1])
dcd_out.open_to_write()

# header
dcd.read_header()
header = dcd.get_header()
dcd_out.set_header(header)
if dcd_out._header.format == 'cafemol':
    dcd_out._header.nmp_real = ID_DOM_END - ID_DOM_INI + 1
    dcd_out._header.nunit_real = 1
    dcd_out._header.lunit2mp = []
    dcd_out._header.lunit2mp.append(ID_DOM_END - ID_DOM_INI + 1)
dcd_out.write_header()

#nmp = header.nmp_real

def unwrap_PBC(d):

    pre = d[0][:]
    add = [0.0] * 3
    for ixyz, xyz in enumerate(d):
        for i in range(3):
            x = xyz[i] + add[i]
            if x - pre[i] > MAXD:
                x      += - BOXSIZE
                add[i] += - BOXSIZE
            elif x - pre[i] < -MAXD:
                x      +=   BOXSIZE
                add[i] +=   BOXSIZE
            pre[i] = x

            d[ixyz][i] = x

    return 


#iframe = 0
while dcd.has_more_data() :
    data = dcd.read_onestep()
    
    ##########################################################
    unwrap_PBC(data[ID_DOM_INI:ID_DOM_END+1])

    #if L[0] > BOXSIZE or L[1] > BOXSIZE or L[2] > BOXSIZE:
    #    print ('Warning: L exceeds BOXSIZE at frame %i' % iframe)

    dcd_out.write_onestep(data[ID_DOM_INI:ID_DOM_END+1]) 

    #iframe += 1
    
dcd.close()
dcd_out.close()
