#!/usr/bin/env python
'''
@author: Naoto Hori
'''

import sys
from lop.file_io.dcd import DcdFile

if (len(sys.argv) != 6):
    print('Usage: % SCRIPT [input dcd] [box x] [box y] [box z] [output dcd')
    print('        e.g. dcd_add_boxinfo.py input.dcd 500.0 500.0 500.0 output.dcd')
    sys.exit(2)
    
box = [float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])]

dcd = DcdFile(sys.argv[1])
dcd.open_to_read()

dcd_out = DcdFile(sys.argv[-1])
dcd_out.open_to_write()

dcd.read_header()
dcd_out.set_header(dcd.get_header())
dcd_out._header.with_unit_cell = True
dcd_out._header.unit_cell_xyz = box
dcd_out._header.nset = 1
dcd_out._header.nstep = 1
dcd_out._header.unit_cell_abc = [0.0, 0.0, 0.0]
dcd_out.write_header()

while dcd.has_more_data():

    dcd_out.write_onestep( dcd.read_onestep() )

dcd.close()
dcd_out.close()
