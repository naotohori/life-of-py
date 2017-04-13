#!/usr/bin/env python

'''
Suppose RNA starts with Sugar.
S, P, B, S, P, B ... so obn
'''

import sys
import numpy as np

if len(sys.argv) != 3:
    print 'Usage: SCRIPT [# mp] [output]'
    sys.exit(2)

nmp = int(sys.argv[1])
f_out = open(sys.argv[-1], 'w')

## mask index starts from 1
mask = np.empty( (nmp+1, nmp+1) )
mask.fill(1)

def mark( target ):
    if 0 < target and target <= nmp:
        mask[ imp, target ] = -1

for imp in range(1, nmp+1):
    res = imp % 3

    if res == 1:  # S
        mark( imp - 1 )
        mark( imp     )
        mark( imp + 1 )
        mark( imp + 2 )

    elif res == 2: # B
        mark( imp - 1 )
        mark( imp     )

    elif res == 0: # P
        mark( imp - 2 )
        mark( imp     )
        mark( imp + 1 )

    else:
        print 'it should not occur....'
        sys.exit(2)

for imp in range(1, nmp+1):
    for jmp in range(1, nmp+1):
        f_out.write(' %2i' % mask[imp, jmp] )
    f_out.write('\n')

