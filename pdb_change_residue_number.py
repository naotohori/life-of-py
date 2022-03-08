#!/usr/bin/env python
'''
@author: Naoto Hori
'''

import sys
from lop.file_io.pdb import PdbFile

if len(sys.argv) != 5 :
    print('')
    print(' Usage: SCRIPT [input PDB] [The residue number of the first a.a.] [output PDB] [log file]')
    print('')
    sys.exit(2)
    
pdb_in = PdbFile(sys.argv[1])
pdb_in.open_to_read()
chains = pdb_in.read_all()
pdb_in.close()

res_id = int(sys.argv[2])
f_log = open(sys.argv[4], 'w')
f_log.write('#original -> new\n')
for c in chains :
    print(sys.argv[1], '#residues', len(c.residues))
    for r in c.residues :
        f_log.write('%i %s -> %i\n' % (r.atoms[0].res_seq, r.atoms[0].ins_code, res_id))
        for a in r.atoms :
            a.res_seq = res_id
            a.ins_code = '    '
        res_id += 1
f_log.close()
        
pdb_out = PdbFile(sys.argv[3])
pdb_out.open_to_write()
pdb_out.write_all(chains)
pdb_out.close()

