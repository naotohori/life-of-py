#!/usr/bin/env python

import sys
from lop.file_io.pdb import PdbFile
from py_gauss_jordan import gauss_jordan
from numpy import float64, zeros

if (len(sys.argv) < 6) or (len(sys.argv) % 2 != 0) :
    print('')
    #print ' Usage: SCRIPT [input pca file] [input pdb file] [(ID begin, ID end) ...] [output pca file]'
    print(' Usage: SCRIPT [input pca file] [input pdb file] [(ID begin, ID end) ...]')
    print('')
    sys.exit(2)

f_in = open(sys.argv[1], 'r')
f_pdb = PdbFile(sys.argv[2])
f_pdb.open_to_read()
chains = f_pdb.read_all()
f_pdb.close()
f_out = open(sys.argv[-1], 'w')

mps = []
for i_pair in range((len(sys.argv) - 3) / 2) :
    id_begin = int(sys.argv[3 + i_pair * 2])
    id_end = int(sys.argv[3 + i_pair * 2 + 1])
    for imp in range(id_begin, id_end + 1) :
        mps.append(imp)
nmp = len(mps)

imp = 0
i_xyz = 0
vecs = []
vec = []
for line in f_in :
    if line.find('#') != -1 :
        continue
    i_xyz += 1
    if i_xyz <= 2: # x,y
        vec.append(float(line.strip()))
    elif i_xyz == 3:
        vec.append(float(line.strip()))
        vecs.append(vec)
        vec = []
        imp += 1
        i_xyz = 0

# coordinates
xyzs = []
imp = 0
for c in chains :
    for imp_chain in range(c.num_atom()) :
        imp += 1
        if imp in mps :
            xyzs.append(c.get_atom(imp_chain).xyz)

#print 'debug: len(mps)  =',len(mps)
#print 'debug: len(vecs) =',len(vecs)
#print 'debug: len(xyzs) =',len(xyzs)
center = [0.0, 0.0, 0.0]
for xyz in xyzs:
    center[0] += xyz.x
    center[1] += xyz.y
    center[2] += xyz.z
print('center=')
print([value / len(xyzs) for value in center])

# calculate translation
translate = [0.0, 0.0, 0.0]
for vec in vecs:
    translate[0] += vec[0]
    translate[1] += vec[1]
    translate[2] += vec[2]
print('translate(BEFORE)=')
print(translate)

# remove translation
translate[0] /= float(nmp)
translate[1] /= float(nmp)
translate[2] /= float(nmp)
for vec in vecs :
    vec[0] -= translate[0]
    vec[1] -= translate[1]
    vec[2] -= translate[2]

# calculate translation
translate = [0.0, 0.0, 0.0]
for vec in vecs:
    translate[0] += vec[0]
    translate[1] += vec[1]
    translate[2] += vec[2]
print('translate(AFTER)=')
print(translate)

# calculate angular momentum
anglmt = [0.0, 0.0, 0.0]
for i in range(len(mps)) :
    anglmt[0] += xyzs[i].y * vecs[i][2] - xyzs[i].z * vecs[i][1]
    anglmt[1] += xyzs[i].z * vecs[i][0] - xyzs[i].x * vecs[i][2]
    anglmt[2] += xyzs[i].x * vecs[i][1] - xyzs[i].y * vecs[i][0]

print('angular momentum(BEFORE)=')
print(anglmt)

txx = 0.0
txy = 0.0
txz = 0.0
tyy = 0.0
tyz = 0.0
tzz = 0.0
for xyz in xyzs :
    txx += xyz.x * xyz.x
    txy += xyz.x * xyz.y
    txz += xyz.x * xyz.z
    tyy += xyz.y * xyz.y
    tyz += xyz.y * xyz.z
    tzz += xyz.z * xyz.z

rot = zeros( (3,3), order='f', dtype=float64)
rot[0,0] = tyy + tzz
rot[1,0] = -txy
rot[2,0] = -txz
rot[0,1] = -txy
rot[1,1] = txx + tzz
rot[2,1] = -tyz
rot[0,2] = -txz
rot[1,2] = -tyz
rot[2,2] = txx + tyy

#print 'rot'
#print rot

etator = zeros((3,3),order='f',dtype=float64)
etator[0,0] = 1.0
etator[1,1] = 1.0
etator[2,2] = 1.0

ier = gauss_jordan(rot, etator, 3)
if ier == -1 :
    print('util_ppgauss is failed')
    sys.exit(2)
    
#print 'etator'
#print etator
    
cx = etator[0,0] * anglmt[0] + etator[0,1] * anglmt[1] + etator[0,2] * anglmt[2]
cy = etator[1,0] * anglmt[0] + etator[1,1] * anglmt[1] + etator[1,2] * anglmt[2]
cz = etator[2,0] * anglmt[0] + etator[2,1] * anglmt[1] + etator[2,2] * anglmt[2]
#print 'cx,cy,cz'
#print cx, cy, cz

for (i,vec) in enumerate(vecs):
    vec[0] += cz * xyzs[i].y -cy * xyzs[i].z
    vec[1] += cx * xyzs[i].z -cz * xyzs[i].x
    vec[2] += cy * xyzs[i].x -cx * xyzs[i].y

# calculate angular momentum
anglmt = [0.0, 0.0, 0.0]
for i in range(len(mps)) :
    anglmt[0] += xyzs[i].y * vecs[i][2] - xyzs[i].z * vecs[i][1]
    anglmt[1] += xyzs[i].z * vecs[i][0] - xyzs[i].x * vecs[i][2]
    anglmt[2] += xyzs[i].x * vecs[i][1] - xyzs[i].y * vecs[i][0]

print('angular momentum(AFTER)=')
print(anglmt)

for vec in vecs:
    f_out.write('%30.20f\n' % vec[0])
    f_out.write('%30.20f\n' % vec[1])
    f_out.write('%30.20f\n' % vec[2])
