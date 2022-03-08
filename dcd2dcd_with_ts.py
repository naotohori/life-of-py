#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2013/08/20
@author: Naoto Hori
'''

import sys
from lop.file_io.dcd import DcdFile
from lop.file_io.ts import TsFile

if len(sys.argv) != 4:
    print('Usage: % SCRIPT [input DCD] [ts file] [output DCD]')
    sys.exit(2)

# Open DCD and read the header
dcd_filename = sys.argv[1]
dcd = DcdFile(dcd_filename)
dcd.open_to_read()
dcd.read_header()

header = dcd.get_header()
#header.istart = header.nstep_save * (frame_begin - 1)
#header.istart = header.istart + header.nstep_save * frame_begin
#header.nset = int(frame_num / frame_stride)
#header.nset = int((frame_end-frame_begin)/frame_stride) + 1
#header.nstep = header.nstep_save * (frame_num - 1)
#header.nstep_save = header.nstep_save * frame_stride

# Open DCD and read the header
dcd_out = DcdFile(sys.argv[-1])
dcd_out.open_to_write()
dcd_out.set_header(header)
dcd_out.write_header()

# Open TS and read the header
ts = TsFile(sys.argv[2])
ts.open_to_read()
ts.read_header()

# read and write
#imodel = 0
#i_org = 0
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
        
        dcd_out.write_onestep(dcd.read_onestep())
        #imodel += 1
    else:
        dcd.skip_onestep()
    
    #i_org += 1

ts.close()
dcd.close()
dcd_out.close()
