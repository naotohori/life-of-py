#!/usr/bin/env python

'''
Created on 2016/06/17
@author: Naoto Hori
'''

import sys
from .file_io.pdb import PdbFile

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: %SCRIPT [input PDB file] [output]')
        sys.exit(2)

f_pdb = PdbFile(sys.argv[1])
f_pdb.open_to_read()
chains = f_pdb.read_all()
f_pdb.close()

if len(chains) != 1:
    print('len(chains) != 1')
    sys.exit(2)

f_out = open(sys.argv[-1],'w')

num_atom = chains[0].num_atom()

for i in range(num_atom):
    xyz_i = chains[0].get_atom(i).xyz
    for j in range(i+1,num_atom):
        xyz_j = chains[0].get_atom(j).xyz
        d = xyz_j.distance(xyz_i)
        f_out.write('%5i %5i %10.2f\n' % (i+1,j+1,d))
