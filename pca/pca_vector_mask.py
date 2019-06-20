#!/usr/bin/env python

import sys

if len(sys.argv) < 5 or len(sys.argv)%2 != 1:
    print ("\n Usage: SCRIPT [input vec file] [(mask ID begin, mask ID end) , ...] [output vec file]\n")
    sys.exit(2)

f_in = open(sys.argv[1], 'r')
f_out = open(sys.argv[-1], 'w')

# input for ID pairs
id_pairs = []
n = 1
for arg in sys.argv[2:-1] :
    if n == 1:
        tp = (int(arg),)
    else :
        id_pairs.append(tp + (int(arg),))
    n *= -1

mask_id = []
for pair in id_pairs:
    mask_id.extend(list(range(pair[0], pair[1]+1))) 

i_xyz = 0
i_mp = 0
for line in f_in :
    if line.find('#') != -1:
        f_out.write(line)
        continue
    i_xyz += 1
    if i_xyz == 1 :
        i_mp += 1
    elif i_xyz == 3:
        i_xyz = 0
        
    if i_mp in mask_id :
        f_out.write('%30.20E\n' % 0.0)
    else: 
        f_out.write(line)
    
f_in.close()
f_out.close()
