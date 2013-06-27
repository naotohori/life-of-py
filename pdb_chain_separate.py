#!/usr/bin/env python

import sys
from file_io.pdb import PdbFile

if len(sys.argv) != 3:
    print ('\n Usage: SCRIPT [input PDB file] [output DIR (with/without prefix)]\n')
    sys.exit(2)
    
f_pdb = PdbFile(sys.argv[1])
f_pdb.open_to_read()
chains = f_pdb.read_all()
f_pdb.close()

print ('%i chains' %len(chains))
for c in chains :
    chain_id = c.get_atom(0).chain_id
    if chain_id == ' ' :
        continue
    filename = sys.argv[2] + chain_id + '.pdb'
    f = open(filename,'w')
    f = PdbFile(filename)
    f.open_to_write()
    f.write_all([c])
