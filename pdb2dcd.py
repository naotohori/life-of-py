#!/usr/bin/env python

from lop.file_io.pdb import PdbFile
from lop.file_io.dcd import DcdFile, DcdHeader

import sys

if len(sys.argv) != 3:
    print('Usage: SCRIPT [pdb] [dcd]')
    sys.exit(2)

pdb = PdbFile(sys.argv[1])
pdb.open_to_read()
chains = pdb.read_all()
pdb.close()

header = DcdHeader()
header.nset = 1
header.istart = 1
header.nstep_save = 1
header.nstep = 1
header.nunit_real = len(chains)
header.delta = 0.0
header.title = [' '*80,' '*80]
header.tempk = 0.0
header.lunit2mp = []
imp = 0
xyzs = []
for c in chains:
    for r in c.residues:
        for a in r.atoms:
            imp += 1
            xyzs.append([a.xyz.x, a.xyz.y, a.xyz.z])
    header.lunit2mp.append(imp)
header.nmp_real = imp

dcd = DcdFile(sys.argv[-1])
dcd.open_to_write()
dcd.set_header(header)
dcd.write_header()
dcd.write_onestep(xyzs)
dcd.close()
