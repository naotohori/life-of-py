#!/usr/bin/env python
#vim:fileencoding=UTF-8

from cafysis.elements.coord import Coord
from math import hypot, atan2, sqrt
from cafysis.mtx_coord_transform import mtx_coord_transform
from scipy.constants.constants import pi

Q = Coord(-2,-2*sqrt(3.0), - 4*sqrt(3.0))
Qxy = hypot(Q.x,Q.y)
print "Qxy",Qxy

#Qと原点を結ぶ線Lと垂直な、xy平面上の軸を表す単位ベクトルR
R = Coord()
R.x =  Q.y / Qxy
R.y = -Q.x / Qxy
R.z = 0.0
print "R",R.get_as_tuple()

phi = atan2(Qxy,Q.z)
print "phi",phi, "(", phi/pi*180.0, "deg)"

rotateQ = mtx_coord_transform()
rotateQ.rotate(R.x, R.y, R.z, phi) 
 
Q.transform(rotateQ.mtx)
print "Q",Q.get_as_tuple()

#ここまでで、QがZ軸上にのった

#NをZ軸まわりに回転させ、Nx>0でNy=0へもっていく。
N = Coord(3,-3,-3)
theta = -atan2(N.y,N.x)
print "theta",theta, "(", theta/pi*180.0, "deg)"
rotateN = mtx_coord_transform()
rotateN.rotate_z(theta)

N.transform(rotateN.mtx)
print "N",N.get_as_tuple()