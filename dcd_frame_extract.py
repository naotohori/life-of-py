#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2011/03/17
Added stride mode 2011/07/29
@author: Naoto Hori
'''

import sys
from lop.file_io.dcd import DcdFile

if (not len(sys.argv) in (5, 6)):
    print (
'''
 Usage: % SCRIPT [input DCD] [beginning (0)] [end] [output DCD]
        % SCRIPT [input DCD] [beginning (0)] [end] [stride] [output DCD]
        
 In this script, frame number starts with "0".
    
  step   DCD frame number
     0        0
   500        1
  1000        2
  1500        3
    ..        .
 10000       20
 
 beginning:0          beginning:0        beginning:2     beginning:1
 end:5                end:20             end:10          end:10
 stride:1(default)    stride:2           stride:2        stride:2
     0        0           0          0         1000      0        500      0
   500        1        1000          1         2000      1       1500      1
  1000        2        2000          2         3000      2       2500      2
  1500        3        3000          3         4000      3       3500      3
  2000        4         ..           .         5000      4       4500      4
  2500        5        9000          9
                      10000         10
 #frame=    6                   11                 5               5

 frame_num   6                  21                 9              10
'''
    )
    sys.exit(2)

frame_begin = int(sys.argv[2])
if frame_begin < 0:
    print('Error: beginning frame should not less than 0')
    sys.exit(2)
frame_end = int(sys.argv[3])
if (len(sys.argv) == 5) :
    frame_stride = 1
else:
    frame_stride = int(sys.argv[4])
    if frame_stride <= 0:
        print('The frame stride is invalid')
        sys.exit(2)
        
frame_num = frame_end - frame_begin + 1
if frame_num < 1 :
    print('The number of frames is invalid.')
    sys.exit(2)

dcd = DcdFile(sys.argv[1])
dcd.open_to_read()
if (len(sys.argv) == 5):
    dcd_out = DcdFile(sys.argv[4])
else:
    dcd_out = DcdFile(sys.argv[5])
dcd_out.open_to_write()

# header
dcd.read_header()
header = dcd.get_header()
##header.istart = header.nstep_save * (frame_begin - 1)
#header.istart = header.istart + header.nstep_save * frame_begin
##header.nset = int(frame_num / frame_stride)
#header.nset = int((frame_end-frame_begin)/frame_stride) + 1
#header.nstep = header.nstep_save * (frame_num - 1)
#header.nstep_save = header.nstep_save * frame_stride
dcd_out.set_header(header)
dcd_out.write_header()

def error_no_data() :
    print('The number of frames is invalid.')
    print('Header information:')
    dcd.show_header()
    sys.exit(2)

# skip
dcd.skip(frame_begin - 1)

# read and write
icount = -1
for i in range(frame_num) :
    if not dcd.has_more_data() :
        error_no_data()
    icount += 1
    if icount % frame_stride == 0 :
        icount = 0
    else :
        dcd.skip(1)
        continue

    data = dcd.read_onestep()

    if dcd_out._header.with_unit_cell:
        dcd_out._header.unit_cell_xyz = dcd._header.unit_cell_xyz
        dcd_out._header.unit_cell_abc = dcd._header.unit_cell_abc

    dcd_out.write_onestep(data)

dcd.close()
dcd_out.close()
