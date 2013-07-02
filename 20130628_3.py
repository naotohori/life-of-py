#!/usr/bin/env python
#vim:fileencoding=UTF-8

'''
RN20:107
20130628_2.pyで、原点にて配向を揃えたDCDに対し、
phi, psi（緯度、経度）を計算する。

いま、chain1=Fus3, chain2=Ste7であり、IDに注意が必要。
chain1 1-353
chain2 354-868
よって、Ste7のコアドメインは、542-824であり、
 Q189 = 542
 N201 = 672

Fus3ドメインのIDは、1-353
'''

from cafysis.file_io.dcd import DcdFile
import sys
from math import atan2, sqrt, acos

if len(sys.argv) != 3:
    print 'Usage: SCRIPT [input DCD] [output file]'
    sys.exit(2)

#RESIDUE_TOTAL = 868
ID_FUS3_DOM_INI = 1 - 1
ID_FUS3_DOM_END = 353 - 1

dcd = DcdFile(sys.argv[1])
dcd.open_to_read()
file_out = open(sys.argv[2], 'w')

# header
dcd.read_header()
nmp = dcd.get_header().nmp_real

#if nmp != RESIDUE_TOTAL:
#    print 'ERROR: nmp != RESIDUE_TOTAL'
#    sys.exit(2)

while dcd.has_more_data() :
    data = dcd.read_onestep()
    
    
    ##########################################################
    #Fus3(1-353)の重心
    com = [0.0] * 3
    n_com = ID_FUS3_DOM_END - ID_FUS3_DOM_INI + 1
    for i in xrange(ID_FUS3_DOM_INI, ID_FUS3_DOM_END+1):
        com[0] += data[i][0]
        com[1] += data[i][1]
        com[2] += data[i][2]
    
    com = [x/float(n_com) for x in com]
    dist = sqrt(com[0]**2 + com[1]**2 + com[2]**2)
    
    theta = acos(com[2]/dist)
    phi = atan2(com[1], com[0])
    
    file_out.write('%f8.3 %f8.3\n' % (theta,phi))
    
    
dcd.close()
file_out.close()