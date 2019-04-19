#!/usr/bin/env python

from numpy import zeros, ones, float64, dot
import sys

if len(sys.argv) != 4:
    print('Usage: % SCRIPT [PDB filename] [matrix filename] [out PDB filename]')
    sys.exit(2)
    
filename_pdb = sys.argv[1]
filename_mat = sys.argv[2]
filename_out_pdb = sys.argv[3]

# Coord1
pdb = PdbFile(filename_pdb)
pdb.open_to_read()
chains = pdb.read_all()
pdb.close()

num_atom = 0
for chain in chains :
    for residue in chain.residues :
        num_atom += len(residue.atoms)
    
file_mat = file(filename_mat,'r')
mat_lines = file_mat.readlines()
if mat_lines[0][0:7] != '#matrix' :
    print('Error: format of matrix file')
    sys.exit()
    
transform = zeros((4,4),dtype=float64)
for i in range(4) :
    transform[i] = mat_lines[i+1].strip().split()
    
file_mat.close()

for chain in chains :
    for residue in chain.residues:
        for atom in residue.atoms :
            target = ones((4), dtype=float64)
            target[0:3] = atom.xyz.get_as_tuple()
            
            target = dot(transform, target.T)
            
            atom.xyz.x = target[0]
            atom.xyz.y = target[1]
            atom.xyz.z = target[2]
    
pdb = PdbFile(filename_out_pdb)
pdb.open_to_write()
pdb.write_all(chains)
pdb.close()

