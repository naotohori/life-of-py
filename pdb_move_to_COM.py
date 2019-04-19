#!/usr/bin/env python

import sys
from file_pdb import PdbFile
from my_element import Coord

if len(sys.argv) != 3:
    print('Usage: % SCRIPT [input pdb] [output pdb]')
    sys.exit(2)
    
f_pdb_in = PdbFile(sys.argv[1])
f_pdb_in.open_to_read()
chains = f_pdb_in.read_all()
f_pdb_in.close()

com = Coord()
total_atom = 0
for c in chains :
    total_atom += c.num_atom()
    for iatom in range(c.num_atom()) :
        com.move(c.get_atom(iatom).xyz)
        
com.x = - com.x / float(total_atom)
com.y = - com.y / float(total_atom)
com.z = - com.z / float(total_atom)

for c in chains :
    for iatom in range(c.num_atom()) :
        c.get_atom(iatom).xyz.move(com)

f_pdb_out = PdbFile(sys.argv[2])
f_pdb_out.open_to_write()
f_pdb_out.write_all(chains)
f_pdb_out.close()
