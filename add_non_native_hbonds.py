#!/usr/bin/env python

import sys
from lop.file_io.pdb import PdbFile

NNHB_NT_SEP = 5

atoms_AU = 'N6 O4 N1 N3 TER'
atoms_GC = 'N1 N3 N2 O2 O6 N4 TER'
atoms_GU = 'N1 O2 O6 N3 TER'

if len(sys.argv) != 4:
    print("Usage: SCRIPT [input cg.pdb (for seq)] [input hb.native.dat] [output hbond.dat]")
    sys.exit(2)

filepath_pdb = sys.argv[1]
filepath_hb = sys.argv[2]
filepath_out = sys.argv[3]


pdb = PdbFile(filepath_pdb)
pdb.open_to_read()
chains = pdb.read_all()

if len(chains) != 1:
    print("Error: len(chains) != 1, not supported")
    sys.exit(2)

seq = ''
for r in chains[0].residues:
    seq += r.atoms[0].res_name.strip()[-1]

nseq = len(seq)
print(seq)
print('length = %i' % nseq)


f_out = open(filepath_out,'w')
 
native_pairs = []
for l in open(filepath_hb):
    lsp = l.split()
    native_pairs.append( (int(lsp[1]), int(lsp[2])) )

    f_out.write(l)


for i in range(2, nseq):
    nt_i = i + 1

    for j in range(i+NNHB_NT_SEP, nseq):
        nt_j = j + 1

        if (nt_i,nt_j) in native_pairs:
            continue

        if seq[i] == 'G' and seq[j] == 'C':
            f_out.write('G-C    %4i %4i  %s\n' % (nt_i, nt_j, atoms_GC))
        if seq[i] == 'C' and seq[j] == 'G':
            f_out.write('G-C    %4i %4i  %s\n' % (nt_j, nt_i, atoms_GC))
        if seq[i] == 'A' and seq[j] == 'U':
            f_out.write('A-U    %4i %4i  %s\n' % (nt_i, nt_j, atoms_AU))
        if seq[i] == 'U' and seq[j] == 'A':
            f_out.write('A-U    %4i %4i  %s\n' % (nt_j, nt_i, atoms_AU))
        if seq[i] == 'G' and seq[j] == 'U':
            f_out.write('G-U    %4i %4i  %s\n' % (nt_i, nt_j, atoms_GU))
        if seq[i] == 'U' and seq[j] == 'G':
            f_out.write('G-U    %4i %4i  %s\n' % (nt_j, nt_i, atoms_GU))

f_out.close()
