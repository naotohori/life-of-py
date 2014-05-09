#!/usr/bin/env python

import sys

if len(sys.argv) < 4 :
    print ('')
    print (' Usage: SCRIPT [vector file] [factor] [output vector file]')
    print ('')
    sys.exit(2)
    
fs_in = open(sys.argv[1], 'r')
factor = float(sys.argv[2])
fs_out = open(sys.argv[-1], 'w')

vec = []
for line in fs_in :
    if line.find('#') != -1 :
        continue
    vec.append(float(line.strip()))
fs_in.close()

for v in vec :
    fs_out.write('%30.20E\n' % (v*factor, ))
fs_out.close()
