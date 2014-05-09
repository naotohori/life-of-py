#!/usr/bin/env python

import sys

if len(sys.argv) != 7:
    print ('\n Usage: SCRIPT [input mpvec file] [input aln file (distance file)] [cutoff] [imp begin] [imp end] [output vec file]\n')
    sys.exit(2)
    
f_vec = open(sys.argv[1], 'r')
f_dst = open(sys.argv[2], 'r')
cutoff = float(sys.argv[3])
imp_begin = int(sys.argv[4])
imp_end = int(sys.argv[5])
f_out = open(sys.argv[-1], 'w')

imp_out_to_in = {}
for line in f_dst:
    if line.find('#') != -1:
        continue
    linesp = line.split()
    imp_in  = int(linesp[1])  # column 
    imp_out = int(linesp[3])  # column 
    distance = float(linesp[4])
    
    if distance > cutoff :
        continue
    
    if imp_out in imp_out_to_in:
        print ('Error: %i is already exist in imp_out_to_in' % imp_out)
        sys.exit(2)
    imp_out_to_in[imp_out] = imp_in
    
mpvec = {}
for line in f_vec :
    if line.find('#') != -1:
        continue
    linesp = line.split()
    imp = int(linesp[1])
    v = (float(linesp[2]), float(linesp[3]), float(linesp[4]))
    mpvec[imp] = v
    
for imp_out in xrange(imp_begin, imp_end+1) :
    if imp_out in imp_out_to_in :
        imp_in = imp_out_to_in[imp_out]
        if not imp_in in mpvec :
            print ('Error: not imp_in in imp_in_to_vec')
            sys.exit(2)
        v = mpvec[imp_in]
    else :
        v = (0.0, 0.0, 0.0)
    f_out.write('%30.20E\n%30.20E\n%30.20E\n' % v)
    
