#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2016/06/13
@author: Naoto Hori
'''

import sys
from .file_io.pdb import PdbFile

if len(sys.argv) != 3:
    print ('\n Usage: SCRIPT [input PDB file] [output distance file]\n')
    sys.exit(2)
    
f_pdb = PdbFile(sys.argv[1])
f_pdb.open_to_read()
chains = f_pdb.read_all()
f_pdb.close()

if len(chains) != 1:
    print('len(chains) != 1')
    sys.exit(2)

f_out = open(sys.argv[2],'w')

n_res = len(chains[0].residues)

for i in range(n_res):
    for j in range(n_res):
        mindist = 9999999.9
        for ai in chains[0].residues[i].atoms:
            for aj in chains[0].residues[j].atoms:
                d = ai.xyz.distance(aj.xyz)
                if d < mindist:
                    mindist = d
        f_out.write('%5i %5i %10.2f\n' % (i+1,j+1,mindist))
    f_out.write('\n')
