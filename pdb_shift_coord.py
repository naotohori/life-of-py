#!/usr/bin/env python
#vim:fileencoding=UTF-8
'''
Created on 2021/12/05
@author: Naoto Hori
'''

from cafysis.mtx_coord_transform import mtx_crd_transform
from cafysis.file_io.pdb import PdbFile
import sys

if len(sys.argv) != 6:
    print('Usage: SCRIPT [input PDB] [x] [y] [z] [output PDB]')
    sys.exit(2)

shift_x = float(sys.argv[2])
shift_y = float(sys.argv[3])
shift_z = float(sys.argv[4])

pdb = PdbFile(sys.argv[1])
pdb.open_to_read()
pdb_out = PdbFile(sys.argv[-1])
pdb_out.open_to_write()

# header
chains = pdb.read_all()

nmp = 0
data = []
for c in chains:
    nmp += c.num_atom()
    for r in c.residues:
        for a in r.atoms:
            data.append([a.xyz.x,a.xyz.y,a.xyz.z])

##########################################################
mtx_move = mtx_crd_transform()
mtx_move.translation(shift_x, shift_y, shift_z)

mtx_move.do_to_data(data)


i = 0
for c in chains:
    for r in c.residues:
        for a in r.atoms:
            a.xyz.x = data[i][0]
            a.xyz.y = data[i][1]
            a.xyz.z = data[i][2]
            i += 1
    
pdb_out.write_all(chains)
pdb_out.close()
