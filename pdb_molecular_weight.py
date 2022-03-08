#!/usr/bin/env python

import sys
from lop.file_io.pdb import PdbFile
from lop.para.mass import ATOM_MASS

if len(sys.argv) != 2:
    print('Usage: SCRIPT (PDB file)')
    sys.exit(2)

pdb = PdbFile(sys.argv[1])
pdb.open_to_read()
chains = pdb.read_all()
pdb.close()

total = 0.0
for ic, c in enumerate(chains):

    weight = 0.0

    for r in c.residues:
        for a in r.atoms:
            weight += ATOM_MASS[a.name[0]]

    print(f'Chain {ic+1}: {weight}')
    total += weight

print(f'Total: {total}')
