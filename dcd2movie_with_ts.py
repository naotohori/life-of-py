#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2013/08/06
@author: Naoto Hori
'''

import sys
from cafysis.file_io.dcd import DcdFile
from cafysis.file_io.ts import TsFile
from cafysis.file_io.pdb import PdbFile

if len(sys.argv) != 5:
    print('Usage: % SCRIPT [input DCD] [ts file] [reference PDB] [output movie]')
    sys.exit(2)

# Read the reference PDB
pdb = PdbFile(sys.argv[-2])
pdb.open_to_read()
chains = pdb.read_all()
pdb.close()

# Output PDB
pdb = PdbFile(sys.argv[-1])
pdb.open_to_write()

# Open DCD and read the header
dcd_filename = sys.argv[1]
dcd = DcdFile(dcd_filename)
dcd.open_to_read()
dcd.read_header()

# Open TS and read the header
ts = TsFile(sys.argv[2])
ts.open_to_read()
ts.read_header()

# read and write
imodel = 0
i_org = 0
while dcd.has_more_data():
    if not ts.has_more_data():
        print('Not enough data in .ts file (1)')
        sys.exit(2)
        
    tsdata, lines = ts.read_onestep()
    
    # skip step=1
    if tsdata[0][ts.head_col.step] == 1:
        if not ts.has_more_data():
            print('Not enough data in .ts file (2)')
            sys.exit(2)
        tsdata, lines = ts.read_onestep()
    
    if tsdata[0][ts.head_col.e_bridge] == '0.0000000':
        struct = dcd.read_onestep()
        iatom = 0
        for c in chains:
            for r in c.residues:
                for a in r.atoms:
                    a.xyz.put_as_list(struct[iatom])
                    iatom += 1
        imodel += 1
        pdb.modelID = imodel
        pdb.set_remark('ORIGINAL_FILE %s' % (dcd_filename,))
        pdb.set_remark('ORIGINAL_FRAME %i' % (i_org,))
        for iunit, x in enumerate(tsdata):
            if iunit == 0:
                tsline = 'TS_all'
            else:
                tsline = 'TS_%i  ' % (iunit,)
            for y in x:
                tsline += ' ' + y
            pdb.set_remark(tsline)
        pdb.write_all(chains)
    else:
        dcd.skip_onestep()
    
    i_org += 1

ts.close()
dcd.close()
pdb.close()
