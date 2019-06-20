#!/usr/bin/env python
#vim:fileencoding=UTF-8

'''
RN20:107
Ste7のPDBについて、位置と配向を揃える
'''

from cafysis.elements.coord import Coord
from cafysis.mtx_coord_transform import mtx_crd_transform
from cafysis.file_io.pdb import PdbFile
import sys
from math import hypot, atan2

if len(sys.argv) != 2:
    print('Usage: SCRIPT [file]')
    sys.exit(2)

pdb = PdbFile(sys.argv[1])
pdb.open_to_read()
chains = pdb.read_all()

#DEBUG from math import sqrt
#DEBUG Q = Coord(-2,-2*sqrt(3.0), - 4*sqrt(3.0))
# QにSte7のQ189:CAの座標を代入
Q = Coord()
for atom in chains[1].residues[188].atoms:
    if atom.name.find("CA") != -1:
        #(Q.x,Q.y,Q.z) = atom.xyz.get_as_tuple()
        Q = atom.xyz
N = Coord()
for atom in chains[1].residues[200].atoms:
    if atom.name.find("CA") != -1:
        N = atom.xyz

# 動かす前        
print("動かす前")
print("Q:",Q.get_as_tuple())
print("N:",N.get_as_tuple())

##########################################################
#Ste7(189-471)の重心(CAのみで計算）が原点に重なるように並進
com = Coord()
n_com = 0
for (i,r) in enumerate(chains[1].residues):
    if (i+1)>=189 and (i+1)<=471:
        for a in r.atoms:
            if a.name.find("CA") != -1:
                com += a.xyz
                n_com += 1

if n_com != 283:
    print("ERROR: n_com = ",n_com)
    sys.exit(2)

com /= n_com
trans_com = mtx_crd_transform()
trans_com.translation(-com.x, -com.y, -com.z)

for c in chains:
    for r in c.residues:
        for a in r.atoms:
            a.xyz.transform(trans_com.mtx)
            
print("重心を原点へ")
print("Q:",Q.get_as_tuple())
print("N:",N.get_as_tuple())


##########################################################
#重心を固定したまま、QをZ軸上へ持っていく。
#Qと原点を結ぶ線Lと垂直な、xy平面上の軸を表す単位ベクトルR
Qxy = hypot(Q.x,Q.y)
#DEBUG print "Qxy",Qxy

R = Coord()
R.x =  Q.y / Qxy
R.y = -Q.x / Qxy
R.z = 0.0
#DEBUG print "R",R.get_as_tuple()

phi = atan2(Qxy,Q.z)
#DEBUG from scipy.constants.constants import pi
#DEBUG print "phi",phi, "(", phi/pi*180.0, "deg)"

rotateQ = mtx_crd_transform()
rotateQ.rotate(R.x, R.y, R.z, phi) 
 
#DEBUG Q.transform(rotateQ.mtx)
#DEBUG print "Q",Q.get_as_tuple()

for c in chains:
    for r in c.residues:
        for a in r.atoms:
            a.xyz.transform(rotateQ.mtx)
#ここまでで、QがZ軸上にのった
print("QをZ軸上へ")
print("Q:",Q.get_as_tuple())
print("N:",N.get_as_tuple())


##########################################################
#NをZ軸まわりに回転させ、Nx>0でNy=0な位置へもっていく。
#DEBUG N = Coord(3,-3,-3)
theta = -atan2(N.y,N.x)
#DEBUG print "theta",theta, "(", theta/pi*180.0, "deg)"
rotateN = mtx_crd_transform()
rotateN.rotate_z(theta)

#DEBUG N.transform(rotateN.mtx)
for c in chains:
    for r in c.residues:
        for a in r.atoms:
            a.xyz.transform(rotateN.mtx)
#DEBUG print "N",N.get_as_tuple()
#ここまでで、NをX軸上へ
print("NをX軸上へ")
print("Q:",Q.get_as_tuple())
print("N:",N.get_as_tuple())

pdb_out = PdbFile("20130628.pdb")
pdb_out.open_to_write()
pdb_out.write_all(chains)
pdb_out.close()