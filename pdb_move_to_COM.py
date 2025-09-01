#!/usr/bin/env python

import sys
from lop.file_io.pdb import PdbFile
from lop.elements.coord import Coord
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--begin", type=int, default=None, help="domain ID begin")
parser.add_argument("--end", type=int, default=None, help="domain ID end")
parser.add_argument('input_pdb', type=PdbFile, help="input PDB")
parser.add_argument('output_pdb', type=PdbFile, help="output PDB")
args = parser.parse_args()

f_pdb_in = args.input_pdb
f_pdb_in.open_to_read()
chains = f_pdb_in.read_all()
f_pdb_in.close()

com = Coord()
natom_total = 0
natom_COM = 0
for c in chains :
    if args.end is not None and natom_total > args.begin:
        break
    for iatom in range(c.num_atom()) :
        natom_total += 1
        if args.begin is not None and natom_total < args.begin:
            continue
        if args.end is not None and natom_total > args.begin:
            continue
        com.move(c.get_atom(iatom).xyz)
        natom_COM += 1

com.x = - com.x / float(natom_COM)
com.y = - com.y / float(natom_COM)
com.z = - com.z / float(natom_COM)

for c in chains :
    for iatom in range(c.num_atom()) :
        c.get_atom(iatom).xyz.move(com)

f_pdb_out = args.output_pdb
f_pdb_out.open_to_write()
f_pdb_out.write_all(chains)
f_pdb_out.close()
