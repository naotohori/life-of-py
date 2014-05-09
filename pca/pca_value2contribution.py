#!/usr/bin/env python

import sys

if len(sys.argv) != 3 :
    print ('')
    print (' Usage: SCRIPT [input pca value file] [output file]')
    print ('')
    sys.exit(2)
    
f_in = open(sys.argv[1],'r')
f_out = open(sys.argv[2], 'w')

values = []
sum = 0.0
for line in f_in :
    if line.find('#') != -1 :
        continue
    value = float(line.strip())
    if value < 0.0 :
        value = 0.0
    values.append(value)
    sum += value
f_in.close()

f_out.write('#contribution\n')
for i,value in enumerate(values):
    f_out.write('%10i %10.5f\n' % (i+1, value/sum*100.0 ,))
    
f_out.write('\n\n')
f_out.write('#contribution (integrated)\n')
integrated = 0.0
for i, value in enumerate(values):
    integrated += value/sum*100.0
    f_out.write('%10i %10.5f\n' % (i+1, integrated))


f_out.close()
