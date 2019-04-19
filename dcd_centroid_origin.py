#!/usr/bin/env python
#vim:fileencoding=UTF-8
'''
Created on 2016/08/05
@author: Naoto Hori

Calculate the centroid for particular part of molecule (defined by ID1 and ID2)
transfer all coordinates so that the centroid comes to the origin (0,0,0).
'''

from cafysis.mtx_coord_transform import mtx_crd_transform
from cafysis.file_io.dcd import DcdFile
import sys

if len(sys.argv) != 5:
    print('Usage: SCRIPT [input DCD] [ID domain begin] [ID domain end] [output DCD]')
    sys.exit(2)

ID_DOM_INI = int(sys.argv[2]) - 1  # 重心を求める際に必要
ID_DOM_END = int(sys.argv[3]) - 1

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

while dcd.has_more_data() :
    data = dcd.read_onestep()
    
    ##########################################################
    #重心(CAのみで計算）が原点に重なるように並進
    com = [0.0] * 3
    n_com = ID_DOM_END - ID_DOM_INI + 1
    for i in range(ID_DOM_INI, ID_DOM_END+1):
        com[0] += data[i][0]
        com[1] += data[i][1]
        com[2] += data[i][2]
    
    com = [-x/float(n_com) for x in com]
    trans_com = mtx_crd_transform()
    trans_com.translation(com[0],com[1],com[2])
    
    #for i in xrange(nmp):
    #    data[i][0:3] = trans_com.do_to_array(data[i])
    trans_com.do_to_data(data)
        
    dcd_out.write_onestep(data) 
    
dcd.close()
dcd_out.close()
