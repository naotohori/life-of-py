#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2011/03/17 (dcd_frame_extract.py)
Added stride mode on 2011/07/29 (dcd_frame_extract.py)
Modified from dcd_frame_extract.py on 2011/12/20
@author: Naoto Hori
'''

import sys
from cafysis.file_io.dcd import DcdFile
from cafysis.file_io.pdb import PdbFile

if (not len(sys.argv) in (6, 7)):
    print 'Usage: % SCRIPT [input DCD] [beginning (1)] [end] [reference PDB] [output movie]'
    print '       % SCRIPT [input DCD] [beginning (1)] [end] [stride] [reference PDB] [output movie]'
    sys.exit(2)

frame_begin = int(sys.argv[2])
frame_end = int(sys.argv[3])
if (len(sys.argv) == 5) :
    frame_stride = 1
else:
    frame_stride = int(sys.argv[4])
    if frame_stride <= 0:
        print 'The frame stride is invalid'
        sys.exit(2)
        
frame_num = frame_end - frame_begin + 1
if frame_num < 1 :
    print 'The number of frames is invalid.'
    sys.exit(2)

dcd = DcdFile(sys.argv[1])
dcd.open_to_read()
#if (len(sys.argv) == 5):
#    dcd_out = DcdFile(sys.argv[4])
#else:
#    dcd_out = DcdFile(sys.argv[5])
#dcd_out.open_to_write()

# header
dcd.read_header()
header = dcd.get_header()
header.istart = frame_begin
header.nset = int(frame_num / frame_stride)
header.nstep = header.nstep_save * (frame_num - 1)
header.nstep_save = header.nstep_save * frame_stride
#dcd_out.set_header(header)
#dcd_out.write_header()

def error_no_data() :
    print 'The number of frames is invalid.'
    print 'Header information:'
    dcd.show_header()
    sys.exit(2)

# skip
dcd.skip(frame_begin - 1)

# read and write
icount = -1
for i in xrange(frame_num) :
    if not dcd.has_more_data() :
        error_no_data()
    icount += 1
    if icount % frame_stride != 0 :
        dcd.skip(1)
        continue
    else :
        icount = 0
#    dcd_out.write_onestep(dcd.read_onestep())

dcd.close()
#dcd_out.close()
