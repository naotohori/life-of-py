#!/usr/bin/env python
#vim:fileencoding=UTF-8

'''
Created on 2016/11/14
@author: Naoto Hori

Ref.
Kevin W Plaxco, Kim T Simons, David Baker
Journal of Molecular Biology
J. Mol. Biol. (1998) 277(4) 985-994
doi:10.1006/jmbi.1998.1645
'''

MIN_IJ = 3
CUTOFF = 8.0

import sys

if len(sys.argv) != 2:
    print('Usage: SCRIPT [pairdist.dat]')
    sys.exit(2)


s = 0
m = 0
n = 0

for l in open(sys.argv[1],'r'):
    lsp = l.split()
    i = int(lsp[0])
    j = int(lsp[1])
    d = float(lsp[2])

    if i > n:
        n = i
    if j > n:
        n = j

    sij = abs(i-j)

    if sij < MIN_IJ:
        continue
    if d >= CUTOFF:
        continue

    s += sij
    m += 1

print('s: ',s)
print('m: ',m)
print('n: ',n)
print('ACO: ', s / float(m))
print('RCO: ', s / float(m*n))
