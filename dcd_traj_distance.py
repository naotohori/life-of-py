#!/usr/bin/env python
'''
Created on 2016/08/05
@author: Naoto Hori
'''

import sys
import math
from lop.file_io.dcd import DcdFile

if len(sys.argv) != 5:
    print(' Usage: % SCRIPT [input DCD] [ID1] [ID2] [output PDB] ')
    sys.exit(2)
    
dcd = DcdFile(sys.argv[1])
id1 = int(sys.argv[2]) - 1
id2 = int(sys.argv[3]) - 1

dcd.open_to_read()
dcd.read_header()

f_out = open(sys.argv[-1],'w')
#nframe = 0
while dcd.has_more_data() :
    data = dcd.read_onestep()
    #nframe += 1
    d = math.sqrt( (data[id1][0] - data[id2][0]) ** 2
                  +(data[id1][1] - data[id2][1]) ** 2
                  +(data[id1][2] - data[id2][2]) ** 2 )
    f_out.write('%.2f\n' % d)
dcd.close()

f_out.close()
