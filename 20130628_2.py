#!/usr/bin/env python
#vim:fileencoding=UTF-8

'''
RN20:107
20130628_1.pyでは、Ste7のPDBについて、位置と配向を揃えた。
このスクリプトでは、まったく同様の手順で、DCDファイルの全構造の位置を揃える。

いま、chain1=Fus3, chain2=Ste7であり、IDに注意が必要。
chain1 1-353
chain2 354-868
よって、Ste7のコアドメインは、542-824であり、
 Q189 = 542
 N201 = 554
 
'''

from cafysis.mtx_coord_transform import mtx_crd_transform
from cafysis.file_io.dcd import DcdFile
import sys
from math import hypot, atan2

if len(sys.argv) != 3:
    print 'Usage: SCRIPT [input DCD] [output DCD]'
    sys.exit(2)

ID_Q = 542 - 1
ID_N = 554 - 1
ID_STE7_DOM_INI = 542 - 1  # 重心を求める際に必要
ID_STE7_DOM_END = 824 - 1
RESIDUE_TOTAL = 868

dcd = DcdFile(sys.argv[1])
dcd.open_to_read()
dcd_out = DcdFile(sys.argv[2])
dcd_out.open_to_write()

# header
dcd.read_header()
header = dcd.get_header()
dcd_out.set_header(header)
dcd_out.write_header()

nmp = header.nmp_real
if nmp != RESIDUE_TOTAL:
    print 'ERROR: nmp != RESIDUE_TOTAL'
    sys.exit(2)

while dcd.has_more_data() :
    data = dcd.read_onestep()
    
    ##########################################################
    #print "動かす前"
    #Q = Coord()
    #Q.put_as_list(data[ID_Q])
    #N = Coord()
    #N.put_as_list(data[ID_N])
    #print "Q:",data[ID_Q]
    #print "N:",data[ID_N]
    
    
    ##########################################################
    #Ste7(189-471)の重心(CAのみで計算）が原点に重なるように並進
    com = [0.0] * 3
    n_com = ID_STE7_DOM_END - ID_STE7_DOM_INI + 1
    #n_com = 0
    for i in xrange(ID_STE7_DOM_INI, ID_STE7_DOM_END+1):
        com[0] += data[i][0]
        com[1] += data[i][1]
        com[2] += data[i][2]
    #    n_com += 1
    
    #if n_com != 283:
    #    print "ERROR: n_com = ",n_com
    #    sys.exit(2)
    
    com = [-x/float(n_com) for x in com]
    trans_com = mtx_crd_transform()
    trans_com.translation(com[0],com[1],com[2])
    
    #for i in xrange(nmp):
    #    data[i][0:3] = trans_com.do_to_array(data[i])
    trans_com.do_to_data(data)
        
    #print "重心を原点へ"
    #Q = Coord()
    #Q.put_as_list(data[ID_Q])
    #N = Coord()
    #N.put_as_list(data[ID_N])
    #print "Q:",data[ID_Q]
    #print "N:",data[ID_N]
    
    
    ##########################################################
    #重心を固定したまま、QをZ軸上へ持っていく。
    #Qと原点を結ぶ線Lと垂直な、xy平面上の軸を表す単位ベクトルR
    
    Q = data[ID_Q][0:3]
    Qxy = hypot(Q[0],Q[1])
    
    R = [0.0] * 3
    R[0] =  Q[1] / Qxy
    R[1] = -Q[0] / Qxy
    
    phi = atan2(Qxy,Q[2])
    
    rotateQ = mtx_crd_transform()
    rotateQ.rotate(R[0], R[1], R[2], phi) 
     
    rotateQ.do_to_data(data)
        
    #print "QをZ軸上へ"
    #Q = Coord()
    #Q.put_as_list(data[ID_Q])
    #N = Coord()
    #N.put_as_list(data[ID_N])
    #print "Q:",data[ID_Q]
    #print "N:",data[ID_N]
    
    
    ##########################################################
    #NをZ軸まわりに回転させ、Nx>0でNy=0な位置へもっていく。
    N = data[ID_N][0:3]
    theta = -atan2(N[1],N[0])
    rotateN = mtx_crd_transform()
    rotateN.rotate_z(theta)
    
    rotateN.do_to_data(data)
        
    #print "NをX軸上へ"
    #Q = Coord()
    #Q.put_as_list(data[ID_Q])
    #N = Coord()
    #N.put_as_list(data[ID_N])
    #print "Q:",data[ID_Q]
    #print "N:",data[ID_N]
    
    dcd_out.write_onestep(data) 
    
dcd.close()
dcd_out.close()