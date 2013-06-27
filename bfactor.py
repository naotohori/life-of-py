#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@author: Naoto Hori
'''
import sys
import math
from numpy import zeros, float64, array
from cafysis.file_io.dcd import DcdFile
from cafysis.file_io.pdb import PdbFile

if len(sys.argv) != 5:
    print 'Usage: % SCRIPT [PDB filename] [DCD filename] [out PDB filename] [out bfactor file]'
    sys.exit(2)
    
filename_pdb = sys.argv[1]
filename_dcd = sys.argv[2]
filename_outpdb = sys.argv[3]
filename_out = sys.argv[4]

# Coord1 fron input PDB
pdb = PdbFile(filename_pdb)
pdb.open_to_read()
chains = pdb.read_all()
pdb.close()

num_atom = 0
for chain in chains :
    num_atom += chain.num_atom()
    
#i = 0
#for chain in ref_chains :
#    for residue in chain.residues:
#        for atom in residue.atoms :
#            (ref[0][i], ref[1][i], ref[2][i]) = atom.xyz.get_as_tuple()
#            i += 1
#
#ref_idx = []
#pre_idx = []

# all to all
#for i in range(num_atom) :
#    ref_idx.append(i + 1)
#    pre_idx.append(i + 1)
    
dcd = DcdFile(filename_dcd)
dcd.open_to_read()
dcd.read_header()

bf = zeros((num_atom,3), dtype=float64)
bf_2 = zeros((num_atom,3), dtype=float64)

k = 0
while dcd.has_more_data() :
    k += 1
    data = array(dcd.read_onestep())
    bf += data
    bf_2 += data**2
    
dcd.close()

bf = bf_2 / k - (bf / k) ** 2
bf_ave = bf.sum() * 8.0*math.pi*math.pi/3.0 / float(num_atom)

bf_out = file(filename_out,'w')
iatom = 0
coef_pdb_bfactor = 10.0
for chain in chains :
    for residue in chain.residues :
        for i in xrange(len(residue.atoms)) :
            bfactor = (bf[iatom,0]+bf[iatom,1]+bf[iatom,2])*8.0*math.pi*math.pi/3.0
            rmsf = math.sqrt(bf[iatom,0]+bf[iatom,1]+bf[iatom,2])
            residue.atoms[i].temp_factor = coef_pdb_bfactor * bfactor / bf_ave
            bf_out.write('%i %f %f %f\n'%(iatom+1,bfactor,bfactor/bf_ave, rmsf))
            iatom += 1
            
# output PDB
pdb_out = PdbFile(filename_outpdb)
pdb_out.open_to_write()
pdb_out.write_all(chains)
pdb_out.close()

bf_out.close()

