#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2011/03/17
Added stride mode 2011/07/29
@author: Naoto Hori
'''

import sys
import argparse

from lop.file_io.dcd import DcdFile

#if (not len(sys.argv) in (5, 6)):
#    print (
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
#    )
#    sys.exit(2)

parser = argparse.ArgumentParser(description='Extract frames from DCD',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--from', dest='frame_begin', default=0,
                    action='store', type=int, help='Beginning frame ID (starting from 0)')

parser.add_argument('--to', dest='frame_end', default=-1,
                    action='store', type=int, help='End frame ID (until the end of file if -1)')

parser.add_argument('--stride', dest='frame_stride', default=1,
                    action='store', type=int, help='Stride of frame')

parser.add_argument('--id-offset', dest='id_offset', default=0,
                    action='store', type=int, help='ID offset for MOdEL')

parser.add_argument('--origin', dest='flg_origin', default=False,
                    action='store_true', help='Move the molecule to the origin')

parser.add_argument('dcd', help='Input DCD file')
parser.add_argument('out', help='Output DCD file')

args = parser.parse_args()


if args.frame_begin < 0:
    print('Error: Beginning frame must be > 0')
    sys.exit(2)

if args.frame_stride <= 0:
    print('Error: Frame stride must be > 0')
    sys.exit(2)

dcd = DcdFile(args.dcd)
dcd.open_to_read()
dcd_out = DcdFile(args.out)
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

num_dcd_frames = dcd.count_frame()

if args.frame_end == -1:
    frame_end = num_dcd_frames - 1

elif args.frame_end > num_dcd_frames - 1:
    print('Warning: There are only %i frames found in the dcd while you specified --to %i' % (num_dcd_frames, args.frame_end))
    print('         Output only %i frames.' % (num_dcd_frames))
    frame_end = num_dcd_frames - 1

else:
    frame_end = args.frame_end

# skip
skipped = dcd.skip_as_many_as_possible_upto(args.frame_begin)
if skipped < args.frame_begin:
    print('Only %i frames could be read from the dcd file, which is less than the one specified as --from %i.' % (skipped, args.frame_begin))
    print('No output file generated.')
    sys.exit(2)

frame_num = frame_end - args.frame_begin + 1
i_orig = args.frame_begin


# read and write
for iframe in range(frame_num) :
    if not dcd.has_more_data() :
        print('The number of frames is invalid.')
        print('Header information:')
        dcd.show_header()
        dcd.close()
        dcd_out.close()
        sys.exit(2)
        
    if iframe % args.frame_stride != 0 :
        dcd.skip(1)
        i_orig += 1
        continue

    struct = dcd.read_onestep()

    if dcd_out._header.with_unit_cell:
        dcd_out._header.unit_cell_xyz = dcd._header.unit_cell_xyz
        dcd_out._header.unit_cell_abc = dcd._header.unit_cell_abc

    """ Move to the origin """
    if args.flg_origin:
        com = [0.0, 0.0, 0.0]
        for v in struct:
            com = [com[i]+v[i] for i in range(3)]

        com = [com[i]/float(len(struct)) for i in range(3)]

        for i in range(len(struct)):
            struct[i] = [struct[i][j] - com[j] for j in range(3)]

    dcd_out.write_onestep(struct)
    i_orig += 1

dcd.close()
dcd_out.close()

