#!/usr/bin/env python

from cafysis.file_io.pdb import PdbFile
from cafysis.elements.pdb import Chain, Residue, Atom
from cafysis.elements.coord import Coord

aa = PdbFile('ru40_aa.pdb')
aa.open_to_read()
chains = aa.read_all()

cg_chains = []

res_id = 0
atom_id = 0
for c in chains:
    for ir, r in enumerate(c.residues):
        xyz_P = Coord()
        xyz_S = Coord()
        xyz_B = Coord()
        for a in r.atoms:
            name = a.name.strip()
            if name[0] == 'H':
                continue
            elif name == 'P':
                xyz_P = a.xyz
            elif name.find("'") != -1:
                xyz_S += a.xyz
            else:
                xyz_B += a.xyz
        res_id += 1
        r_cg = Residue()
        if ir != 0:
            atom_id += 1
            a = Atom()
            a.xyz = xyz_P
            a.res_id = res_id
            a.name = ' P  '
            a.res_name = r.name
        a = Atom()
        a.xyz = xyz_P






