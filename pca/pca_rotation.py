#!/usr/bin/env python

import sys
from file_pdb import PdbFile

if (len(sys.argv) < 5) or (len(sys.argv) % 2 != 1) :
    print('')
    #print ' Usage: SCRIPT [input pca file] [input pdb file] [(ID begin, ID end) ...] [output pull file]'
    print(' Usage: SCRIPT [input pca file] [input pdb file] [(ID begin, ID end) ...]')
    print('')
    sys.exit(2)
    
f_in = open(sys.argv[1], 'r')
f_pdb = PdbFile(sys.argv[2])
f_pdb.open_to_read()
chains = f_pdb.read_all()
f_pdb.close()
#f_out = open(sys.argv[-1], 'w')

mps = []
for i_pair in range((len(sys.argv)-3) / 2) :
    id_begin = int(sys.argv[3+i_pair*2])
    id_end = int(sys.argv[3+i_pair*2+1])
    for imp in range(id_begin, id_end+1) :
        mps.append(imp)
        
imp = 0
i_xyz = 0 
vecs = []
vec = ()
for line in f_in :
    if line.find('#') != -1 :
        continue
    i_xyz += 1
    if i_xyz <= 2: # x,y
        vec += (float(line.strip()) ,)
    elif i_xyz == 3:
        vec += (float(line.strip()) ,)
        vecs.append(vec)
        vec = ()
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

#center = [0.0,0.0,0.0]
#for xyz in xyzs:
#    center[0] += xyz.x
#    center[1] += xyz.y
#    center[2] += xyz.z
#print 'center='
#print [value / len(xyzs) for value in center]

translate = [0.0,0.0,0.0]
for vec in vecs:
    translate[0] += vec[0]
    translate[1] += vec[1]
    translate[2] += vec[2]
print('translation=')
#print [value / len(xyzs) for value in translate]
print(translate)
    
    
anglmt = [0.0, 0.0, 0.0]
for i in range(len(mps)) :
    anglmt[0] += xyzs[i].y * vecs[i][2] - xyzs[i].z * vecs[i][1]
    anglmt[1] += xyzs[i].z * vecs[i][0] - xyzs[i].x * vecs[i][2]
    anglmt[2] += xyzs[i].x * vecs[i][1] - xyzs[i].y * vecs[i][0]
    
print('rotation=')
print(anglmt)
    
    
