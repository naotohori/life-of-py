#!/usr/bin/env python

import sys
from file_pdb import PdbFile
from numpy import float64, zeros

if (len(sys.argv) < 6) or (len(sys.argv) % 2 != 0) :
    print('')
    print(' Usage: SCRIPT [input pca file] [input pdb file] [(ID begin, ID end) ...] [output pca file]')
    #print ' Usage: SCRIPT [input pca file] [input pdb file] [(ID begin, ID end) ...]'
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

for vec in vecs:
    f_out.write('%30.20f\n' % vec[0])
    f_out.write('%30.20f\n' % vec[1])
    f_out.write('%30.20f\n' % vec[2])