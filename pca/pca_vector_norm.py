#!/usr/bin/env python

import sys
import math
from numpy.linalg import norm

if not len(sys.argv) in (2,3):
    print ('\n Usage: SCRIPT [input vec file] (just to confirm)')
    print ('        SCRIPT [input vec file] [output vec file] (normalization)\n')
    sys.exit(2)
    
f_vec = open(sys.argv[1], 'r')
vec = []
for line in f_vec :
    if line.find('#') != -1 :
        continue
    vec.append(float(line))
f_vec.close()
    
print(('dimension = %i' % len(vec)))

#norm = 0.0
#for v in vec :
#    norm += v * v
#    
#norm = math.sqrt(norm)

for i in range(len(vec)/3) :
    a = vec[3*i]**2 + vec[3*i+1]**2 + vec[3*i+2]**2
    print(('%i %30.20f' % (3*i, math.sqrt(a))))
    
value_norm = norm(vec)
print(('norm = %f' % value_norm))

if len(sys.argv) == 3:
    f_out = open(sys.argv[2], 'w')
    for v in vec:
        f_out.write('%30.20E\n' % (v/value_norm ,))
    f_out.close()
    