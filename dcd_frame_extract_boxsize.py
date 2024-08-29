#!/usr/bin/env python
# Huong T Vu Aug 28th 2024

import os
import subprocess
import sys
import math
from lop import dcd_frame_count as dcd_count
from lop.file_io.dcd import DcdFile
from lop.file_io.pdb import PdbFile

if __name__ == "__main__":

    if len(sys.argv) != 7:
        print ('Usage: SCRIPT [dcd file] [boxsize_x 300.000] [boxsize_y 300.000] [boxsize_z 300.000] [pdb reference] [out prefix]')
        sys.exit(2)

dcd = DcdFile(sys.argv[1])
dcd.open_to_read()
dcd.read_header()
num_dcd_frames = dcd.count_frame()

bx = round(float(sys.argv[2]),3)
by = round(float(sys.argv[3]),3)
bz = round(float(sys.argv[4]),3)

# Read the reference PDB
pdb_name=sys.argv[5]
pdb = PdbFile(pdb_name)
pdb.open_to_read()
chains = pdb.read_all()
pdb.close()

frame_l = [] 
i = 0
while dcd.has_more_data():
    struct = dcd.read_onestep()
    if dcd._header.with_unit_cell:
        boxsize = dcd._header.unit_cell_xyz
        if bx == round(boxsize[0],3) and by == round(boxsize[1],3) and bz == round(boxsize[2],3):
            frame_l.append(i)
            iatom = 0
            for c in chains:
                for r in c.residues:
                    for a in r.atoms:
                        a.xyz.put_as_list(struct[iatom])
                        iatom += 1
            pdbout = PdbFile(sys.argv[6]+'_'+str(i)+'.pdb')
            pdbout.open_to_write()
            pdbout.write_all(chains)
            pdbout.close()
    else:
        print('No box info.')
    i += 1

if len(frame_l) == 0:
    print('found no such box size.')
else:
    print('found '+str(len(frame_l))+' frames. Write out as pdb files.')

