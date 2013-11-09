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

if len(sys.argv) != 7:
    print 'Usage: SCRIPT [input DCD] [ID domain begin] [ID domain end] [ID for North] [ID for Greenwich] [output DCD]'
    sys.exit(2)

ID_DOM_INI = int(sys.argv[2]) - 1  # 重心を求める際に必要
ID_DOM_END = int(sys.argv[3]) - 1
ID_N = int(sys.argv[4]) - 1  # North
ID_G = int(sys.argv[5]) - 1  # Greenwich

#ID_N = 542 - 1
#ID_G = 554 - 1
#ID_DOM_INI = 542 - 1  # 重心を求める際に必要
#ID_DOM_END = 824 - 1

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
    #print "動かす前"
    #N = Coord()
    #N.put_as_list(data[ID_N])
    #G = Coord()
    #G.put_as_list(data[ID_G])
    #print "N:",data[ID_N]
    #print "G:",data[ID_G]
    
    
    ##########################################################
    #重心(CAのみで計算）が原点に重なるように並進
    com = [0.0] * 3
    n_com = ID_DOM_END - ID_DOM_INI + 1
    for i in xrange(ID_DOM_INI, ID_DOM_END+1):
        com[0] += data[i][0]
        com[1] += data[i][1]
        com[2] += data[i][2]
    
    com = [-x/float(n_com) for x in com]
    trans_com = mtx_crd_transform()
    trans_com.translation(com[0],com[1],com[2])
    
    #for i in xrange(nmp):
    #    data[i][0:3] = trans_com.do_to_array(data[i])
    trans_com.do_to_data(data)
        
    #print "重心を原点へ"
    #N = Coord()
    #N.put_as_list(data[ID_N])
    #G = Coord()
    #G.put_as_list(data[ID_G])
    #print "N:",data[ID_N]
    #print "G:",data[ID_G]
    
    
    ##########################################################
    #重心を固定したまま、NをZ軸上へ持っていく。
    #Nと原点を結ぶ線Lと垂直な、xy平面上の軸を表す単位ベクトルR
    
    N = data[ID_N][0:3]
    Nxy = hypot(N[0],N[1])
    
    R = [0.0] * 3
    R[0] =  N[1] / Nxy
    R[1] = -N[0] / Nxy
    
    phi = atan2(Nxy,N[2])
    
    rotateN = mtx_crd_transform()
    rotateN.rotate(R[0], R[1], R[2], phi) 
     
    rotateN.do_to_data(data)
        
    #print "NをZ軸上へ"
    #N = Coord()
    #N.put_as_list(data[ID_N])
    #G = Coord()
    #G.put_as_list(data[ID_G])
    #print "N:",data[ID_N]
    #print "G:",data[ID_G]
    
    
    ##########################################################
    #GをZ軸まわりに回転させ、x>0でy=0な位置へもっていく。
    G = data[ID_G][0:3]
    theta = -atan2(G[1],G[0])
    rotateG = mtx_crd_transform()
    rotateG.rotate_z(theta)
    
    rotateG.do_to_data(data)
        
    #print "GをX軸上へ"
    #N = Coord()
    #N.put_as_list(data[ID_N])
    #G = Coord()
    #G.put_as_list(data[ID_G])
    #print "N:",data[ID_N]
    #print "G:",data[ID_G]
    
    dcd_out.write_onestep(data) 
    
dcd.close()
dcd_out.close()