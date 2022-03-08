#!/usr/bin/env python
'''
Created on 2011/05/25
@author: Naoto Hori
'''

import sys
from lop.file_io.dcd import DcdFile
from lop.file_io.pdb import PdbFile

if len(sys.argv) != 4:
    print(' Usage: % SCRIPT [input DCD] [input PDB] [output PDB] ')
    sys.exit(2)
    
dcd = DcdFile(sys.argv[1])
dcd.open_to_read()
dcd.read_header()
nmp = dcd.get_header().nmp_real

f_pdb = PdbFile(sys.argv[2])
f_pdb.open_to_read()
chains = f_pdb.read_all()
f_pdb.close()

f_out = PdbFile(sys.argv[3])
f_out.open_to_write()

ave = []
for i in range(nmp) :
    xyz = [0.0, 0.0, 0.0]
    ave.append(xyz)
    
nframe = 0
while dcd.has_more_data() :
    data = dcd.read_onestep()
    nframe += 1
    print(nframe)
    for i in range(nmp) :
        ave[i][0] += data[i][0]
        ave[i][1] += data[i][1]
        ave[i][2] += data[i][2]
dcd.close()
    
for xyz in ave :
    xyz[0] /= float(nframe)
    xyz[1] /= float(nframe)
    xyz[2] /= float(nframe)
    
imp = 0
for c in chains:
    for i in range(c.num_atom()) :
        imp += 1
        c.get_atom(i).xyz.x = ave[i][0]
        c.get_atom(i).xyz.y = ave[i][1]
        c.get_atom(i).xyz.z = ave[i][2]

f_out.write_all(chains)
f_out.close()
