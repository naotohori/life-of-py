#!/usr/bin/env python

from lop.file_io.pdb import PdbFile
from lop.elements.pdb import Chain, Residue, Atom
import copy
import sys

Nx = 3
Ny = 3
Nz = 3
Lx = 100.0
Ly = 100.0
Lz = 100.0

pdb = PdbFile('./miR21inh.cg.pdb')
pdb.open_to_read()
chains = pdb.read_all()
pdb.close()

newchains = []

for i in range(Nx):
    for j in range(Ny):
        for k in range(Nz):

            dx = i * Lx
            dy = j * Ly
            dz = k * Lz

            for c in chains:
                c_new = Chain()

                for r in c.residues:
                    r_new = Residue()

                    for a in r.atoms:
                        a_new = copy.deepcopy(a)

                        a_new.xyz.x = a.xyz.x + dx
                        a_new.xyz.y = a.xyz.y + dy
                        a_new.xyz.z = a.xyz.z + dz

                        r_new.push_atom(a_new)

                    c_new.push_residue(r_new)

                newchains.append(c_new)

pdb = PdbFile('miR21inh.cg.27L100.pdb')
pdb.open_to_write()
pdb.write_all(newchains)
pdb.close()
