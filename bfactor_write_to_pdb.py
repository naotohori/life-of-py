#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
@author: Naoto Hori
'''
import sys
from lop.file_io.pdb import PdbFile

if len(sys.argv) != 6:
    print ('\n Usage: SCRIPT [bfactor file] [PDB file] [scale] [Upper] [output PDB file]\n')
    sys.exit(2)

f_bf_in = open(sys.argv[1], 'r')
f_pdb = PdbFile(sys.argv[2])
f_pdb.open_to_read()
chains = f_pdb.read_all()
f_pdb.close()

f_pdb_out = PdbFile(sys.argv[-1])
f_pdb_out.open_to_write()
f_pdb_out._file.write('RECORD # SCRIPT: bfactor_write_to_pdb.py\n')
f_pdb_out._file.write('RECORD # argv[1]: '+sys.argv[1]+'\n')
f_pdb_out._file.write('RECORD # argv[2]: '+sys.argv[2]+'\n')
f_pdb_out._file.write('RECORD # argv[3]: '+sys.argv[3]+'\n')
f_pdb_out._file.write('RECORD # argv[4]: '+sys.argv[4]+'\n')
f_pdb_out._file.write('RECORD # argv[5]: '+sys.argv[5]+'\n')

scale = float(sys.argv[3])
upper = float(sys.argv[4])

bf = []
for line in f_bf_in :
    if line.find('#') != -1 :
        continue
    bf.append(float(line.split()[1]))
    
imp = 0
for c in chains:
    for r in c.residues:
        for a in r.atoms:
            x = bf[imp] * scale
            if x > upper :
                x = upper
            a.temp_factor = x
            imp += 1
            
f_pdb_out.write_all(chains)
f_pdb_out.close()
