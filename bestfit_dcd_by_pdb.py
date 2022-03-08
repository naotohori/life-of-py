#!/usr/bin/env python
'''
@author: Naoto Hori
'''

from lop.file_io.dcd import DcdFile
from lop.file_io.pdb import PdbFile
from Superimpose import superimpose
from numpy import zeros, float64
import sys

if len(sys.argv) != 4:
    print('Usage: % SCRIPT [PDB filename] [DCD filename] [out DCD filename]')
    sys.exit(2)
    
filename_pdb = sys.argv[1]
filename_dcd = sys.argv[2]
filename_out_dcd = sys.argv[3]

# Coord1
pdb = PdbFile(filename_pdb)
pdb.open_to_read()
ref_chains = pdb.read_all()
pdb.close()

num_atom = 0
for chain in ref_chains :
    num_atom += chain.num_atom()
    
ref = zeros((3, num_atom), dtype=float64, order='F')

i = 0
for chain in ref_chains :
    for residue in chain.residues:
        for atom in residue.atoms :
            (ref[0][i], ref[1][i], ref[2][i]) = atom.xyz.get_as_tuple()
            i += 1

dcd = DcdFile(filename_dcd)
dcd.open_to_read()
dcd.read_header()

out_dcd = DcdFile(filename_out_dcd)
out_dcd.open_to_write()
out_dcd.set_header(dcd.get_header())
out_dcd.write_header()

k = 0
while dcd.has_more_data() :
    k += 1
    data = dcd.read_onestep_np()
    rmsd = superimpose(ref, data.T) 

    print(k, rmsd)
    out_dcd.write_onestep(data)
        
dcd.close()
out_dcd.close()

