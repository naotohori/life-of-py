#!/usr/bin/env python

import sys
from file_pdb import PdbFile

if not len(sys.argv) in (4,5) :
    print ('\n make a vector which represent movement from PDB1 to PDB2')
    print ('\n Usage: [PDB 1] [PDB 2] [output mpvec file]')
    print (' Usage: [PDB 1] [PDB 2] [output mpvec file] [output vec file]\n')
    sys.exit(2)
    
f_pdb = PdbFile(sys.argv[1])
f_pdb.open_to_read()
chains1 = f_pdb.read_all()
f_pdb.close()
f_pdb = PdbFile(sys.argv[2])
f_pdb.open_to_read()
chains2 = f_pdb.read_all()
f_pdb.close()

f_out_mpvec = open(sys.argv[3], 'w')
if len(sys.argv) == 5:
    f_out_vec   = open(sys.argv[4], 'w')

if len(chains1) != len(chains2) :
    print ('Error: len(chains1) != len(chains2)')
    sys.exit(2)
    
for ic, c1 in enumerate(chains1):
    c2 = chains2[ic]
    
    if c1.num_atom() != c2.num_atom() :
        print ('Error: c1.num_atom() != c2.num_atom()')
        sys.exit(2)
        
    for i in range(c1.num_atom()):
        a1 = c1.get_atom(i)
        a2 = c2.get_atom(i)
        v = (a2.xyz.x - a1.xyz.x, a2.xyz.y - a1.xyz.y, a2.xyz.z - a1.xyz.z)
        f_out_mpvec.write('%5i %6i %30.20E %30.20E %30.20E\n' % ((ic+1, a1.serial) + v))
        if len(sys.argv) == 5:
            f_out_vec.write('%30.20E\n%30.20E\n%30.20E\n' % v)
    