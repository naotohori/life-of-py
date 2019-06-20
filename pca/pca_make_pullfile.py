#!/usr/bin/env python

import sys

if (len(sys.argv) < 5) or (len(sys.argv) % 2 != 1) :
    print('')
    print(' Usage: SCRIPT [input pca file] [(ID begin, ID end) ...] [output pull file]')
    print('')
    sys.exit(2)
    
f_in = open(sys.argv[1], 'r')
f_out = open(sys.argv[-1], 'w')

mps = []
for i_pair in range((len(sys.argv)-3) / 2) :
    id_begin = int(sys.argv[2+i_pair*2])
    id_end = int(sys.argv[2+i_pair*2+1])
    for imp in range(id_begin, id_end+1) :
        mps.append(imp)
    
imp = 0
i_xyz = 0 
vec = ()
for line in f_in :
    if line.find('#') != -1 :
        continue
    i_xyz += 1
    if i_xyz <= 2: # x,y
        vec += (float(line.strip()) ,)
    elif i_xyz == 3:
        vec += (float(line.strip()) ,)
        f_out.write('%i %f %f %f\n' % ((mps[imp],)+vec) )
        vec = ()
        imp += 1
        i_xyz = 0
    
print(('# of mass points : %i' % imp))