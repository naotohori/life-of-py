#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2011/03/17
Added stride mode 2011/07/29
Large refactoring 2025/08/31
@author: Naoto Hori
'''

import sys
import argparse

from lop.file_io.dcd import DcdFile

'''
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

parser = argparse.ArgumentParser(description='Extract frames from DCD',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--from', dest='frame_begin', default=0,
                    action='store', type=int, help='Beginning frame ID (starting from 0)')

parser.add_argument('--to', dest='frame_end', default=-1,
                    action='store', type=int, help='End frame ID (until EOF if -1)')

parser.add_argument('--stride', dest='frame_stride', default=1,
                    action='store', type=int, help='Stride of frame')

parser.add_argument('--origin', dest='flg_origin', default=False,
                    action='store_true', help='Move the molecule to the origin')

parser.add_argument('dcd', help='Input DCD file')
parser.add_argument('out', help='Output DCD file')

args = parser.parse_args()

if args.frame_begin < 0:
    raise SystemExit('Error: --from must be >= 0')

if args.frame_stride <= 0:
    raise SystemExit('Error: --stride must be > 0')

dcd = DcdFile(args.dcd)
dcd.open_to_read()
dcd_out = DcdFile(args.out)
dcd_out.open_to_write()

# header
dcd.read_header()
header = dcd.get_header()
dcd_out.set_header(header)
dcd_out.write_header()

# Go to the "from" frame
if args.frame_begin:
    dcd.go_to_frame(args.frame_begin)

# Maximum frame (--to)
max_keep = None if args.frame_end < 0 else ((args.frame_end - args.frame_begin) // args.frame_stride + 1)

skip = args.frame_stride - 1

wrote = 0
while True:
    if max_keep is not None and wrote >= max_keep:
        break

    ok = dcd.copy_one_frame(dcd_out, origin=args.flg_origin)
    if not ok:
        break  # EOF
    wrote += 1

    # Skip (stride-1)
    if skip > 0:
        dcd.skip(skip)

dcd.close()
dcd_out.flush()
dcd_out.close()

