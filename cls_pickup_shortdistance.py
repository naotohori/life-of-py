#!/usr/bin/env python

import sys

if len(sys.argv) != 5:
    print ('\n Usage: SCRIPT [distance file] [mpvec file] [cutoff] [output vec file]\n')
    sys.exit(2)
    
f_dist = open(sys.argv[1], 'r')
f_mpvec = open(sys.argv[2],'r')
cutoff = float(sys.argv[3])
f_out = open(sys.argv[-1],'w')

imps_output = []
for line in f_dist :
    linesp = line.split()
    imp = int(linesp[1]) # column
    dist = float(linesp[4])
    if dist <= cutoff :
        imps_output.append(imp)
f_dist.close()

for line in f_mpvec:
    linesp = line.split()
    imp = int(linesp[1])   # column
    if imp in imps_output :
        v = (float(linesp[2]), float(linesp[3]), float(linesp[4]))
    else :
        v = (0.0, 0.0, 0.0)
    f_out.write('%30.20E\n%30.20E\n%30.20E\n' % v)
f_mpvec.close()

f_out.close
