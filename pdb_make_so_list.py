#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2014/07/17
@author: Naoto Hori
'''

import sys
from file_io.pdb import PdbFile

if len(sys.argv) != 3:
    print ('\n Usage: SCRIPT [input PDB file] [output SOF list file]\n')
    sys.exit(2)
    
f_pdb = PdbFile(sys.argv[1])
f_pdb.open_to_read()
chains = f_pdb.read_all()
f_pdb.close()

if len(chains) != 1:
    print 'len(chains) != 1'
    sys.exit(2)

f_out = open(sys.argv[2],'w')

for c in chains:
    for i in range(c.num_atom()):
        xyz_i = c.get_atom(i).xyz
        # i is sugar
        if i % 3 == 0:
            j_start = i+3
        # i is base
        elif i % 3 == 1:
            j_start = i+1
        # i is phosphate
        else:
            j_start = i+2

        for j in range(j_start, c.num_atom()):
            d = c.get_atom(j).xyz.distance(xyz_i)
            f_out.write('%4i %4i %6.2f\n' % (i+1,j+1,d))

f_out.close()
