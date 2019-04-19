#!/usr/bin/env python

import sys

if len(sys.argv) < 4 :
    print ('')
    print (' Usage: SCRIPT [vector file 1] [vector file2] [[file 3] [file 4]...] [output vector file]')
    print ('')
    sys.exit(2)
    
fs_in = []
fs_out = open(sys.argv[-1], 'w')
for arg in sys.argv[1:-1] :
    fs_in.append(open(arg,'r'))
    fs_out.write('# add ' + arg + '\n')

vec = []
for line in fs_in[0] :
    if line.find('#') != -1 :
        continue
    vec.append(float(line.strip()))

n = len(vec)

for f in fs_in[1:] :
    i = 0 
    for line in f :
        if line.find('#') != -1:
            continue
        vec[i] += float(line.strip())
        i += 1
    if i != n :
        print ('')
        print(('Error: i != n; i=%i, n=%i', (i,n)))
        print ('')
        sys.exit(2)

for v in vec :
    fs_out.write('%30.20E\n' % v)
    
    