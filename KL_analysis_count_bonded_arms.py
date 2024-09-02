#!/usr/bin/env python

import sys

NMOL = 16
NARM = 4

THRESHOLD = 4

if len(sys.argv) != 3:
    print('SCRIPT [input] [output]')
    print(' input should be output from sisbp_KL_nbp_in_arms_16X.py')
    sys.exit(2)

fout = open(sys.argv[2], 'w')

for l in open(sys.argv[1]):
    lsp = list(map(int, l.strip().split()))
    
    counter = {0:0, 1:0, 2:0, 3:0, 4:0}

    i = -1
    for imol in range(NMOL):

        # Number of arms that have equal or more than THRESHOLD
        # base pairs (min=0, max=4)
        n = 0 

        for iarm in range(NARM):
            i += 1
            if lsp[i] >= THRESHOLD:
                n += 1
        counter[n] += 1

    for i in range(NARM+1):
        fout.write(f' {counter[i]}')
    fout.write('\n')


