#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2014/07/17
@author: Naoto Hori
'''
import sys
import math
from cafysis.file_io.dcd import DcdFile

if len(sys.argv) != 5:
    print ('\n Usage: SCRIPT [input DCD file] [input SO list] [tolerance distance] [output SO file]\n')
    sys.exit(2)
    
#TOLERANCE = 2.0
TOLERANCE = float(sys.argv[3])

solist = []
for l in  open(sys.argv[2],'r'):
    l = l.split()
    i = int(l[0])
    j = int(l[1])
    d = float(l[2])
    solist.append((i,j,d))

ncon = len(solist)

dcd = DcdFile(sys.argv[1])
dcd.open_to_read()
dcd.read_header()

nmp = dcd.get_header().nmp_real

if (nmp-1)*(nmp-2)/2 != ncon:
    print('Warning: (nmp-1)*(nmp-2)/2 != ncon,  ncon=%i, nmp=%i' % (ncon,nmp))
    #sys.exit(2)

f_out = open(sys.argv[-1],'w')
f_out.write('# tolerance distance = %f\n' % (TOLERANCE,))

while dcd.has_more_data():
    data = dcd.read_onestep()

    n = 0
    for (i,j,d) in solist:
        xyz_i = data[i-1]
        xyz_j = data[j-1]
        dist = math.sqrt ((xyz_i[0] - xyz_j[0]) ** 2
                        + (xyz_i[1] - xyz_j[1]) ** 2
                        + (xyz_i[2] - xyz_j[2]) ** 2)
        if abs(dist - d) <= TOLERANCE:
            n += 1
    ki = n / float(ncon)

    f_out.write('%f\n' % ki)

dcd.close()
f_out.close()
