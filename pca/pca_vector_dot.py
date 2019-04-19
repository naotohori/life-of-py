#!/usr/bin/env python

import sys
from numpy import dot

if not len(sys.argv) in (3,):
    print ('\n Usage: SCRIPT [input vec filer1] [input vec file2]\n')
    sys.exit(2)
    
f_vec1 = open(sys.argv[1], 'r')
f_vec2 = open(sys.argv[2], 'r')

vec1 = []
for line in f_vec1 :
    if line.find('#') != -1 :
        continue
    vec1.append(float(line))
f_vec1.close()

vec2 = []
for line in f_vec2 :
    if line.find('#') != -1 :
        continue
    vec2.append(float(line))
f_vec2.close()

if len(vec1) != len(vec2) :
    print ('Error: two vectors have different dimension.')
    print(('       len(vec1) = %i, len(vec2) =%i\n' %(len(vec1), len(vec2))))
    sys.exit(2)
    
print(('dimension = %i' % len(vec1)))

print(dot(vec1, vec2))
