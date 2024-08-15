#!/usr/bin/env python
'''
Created on 2016/08/05
@author: Naoto Hori
'''

import sys
import math
from lop.file_io.dcd import DcdFile
from lop.torsion import torsion

if len(sys.argv) not in (3, 5):
    print(' Usage: % SCRIPT [input DCD] [ID1] [ID2] [ID3] [ID4] [output] ')
    sys.exit(2)
    
dcd = DcdFile(sys.argv[1])

flg_all = False
if len(sys.argv) == 3:
    flg_all = True

if not flg_all:
    id1 = int(sys.argv[2]) - 1
    id2 = int(sys.argv[3]) - 1
    id3 = int(sys.argv[4]) - 1
    id4 = int(sys.argv[5]) - 1

dcd.open_to_read()
dcd.read_header()

f_out = open(sys.argv[-1],'w')
#nframe = 0
while dcd.has_more_data() :
    data = dcd.read_onestep_np()
    #nframe += 1

    if flg_all:
        nmp = len(data)
    
        for i in range(nmp-3):
            dih = torsion(data[i], data[i+1], data[i+2], data[i+3])
            f_out.write(' %.2f' % dih)
        f_out.write('\n')

    else:
        dih = torsion(data[id1], data[id2], data[id3], data[id4])
        f_out.write('%.2f\n' % dih)


dcd.close()

f_out.close()
