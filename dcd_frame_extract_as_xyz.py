#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2011/03/17 (dcd_frame_extract.py)
Added stride mode on 2011/07/29 (dcd_frame_extract.py)
Modified from dcd_frame_extract_as_movie.py on 2023/01/22
Added flg_4500, 2025/3/16
@author: Naoto Hori
'''

import sys
import argparse

from lop.file_io.dcd import DcdFile

parser = argparse.ArgumentParser(description='Convert DCD to Movie (sequence of pdb format)',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--frame', dest='frame', default=None,
                    action='store', type=int, help='Frame ID (starting from 0, default: last frame)')

parser.add_argument('--id-offset', dest='id_offset', default=0,
                    action='store', type=int, help='ID offset for MODEL')


group = parser.add_mutually_exclusive_group()
group.add_argument('--origin', dest='flg_origin', default=False,
                    action='store_true', help='Move the molecule to the origin')

group.add_argument('--4500', dest='flg_4500', default=False,
                    action='store_true', help='Move the molecule centre to (4500, 4500, 4500)')

parser.add_argument('dcd', help='Input DCD file')
parser.add_argument('xyz', help='Reference xyz file')
parser.add_argument('out', help='Output xyz file')

args = parser.parse_args()


if args.frame is not None and  args.frame < 0:
    print('Error: Frame ID must be >= 0')
    sys.exit(2)

# Read the reference XYZ
seq = []
il = 0
for l in open(args.xyz):
    if il > 1:
        lsp = l.split()
        #xyz.append([float(lsp[1]), float(lsp[2]), float(lsp[3])])
        seq.append(lsp[0])

    if il == 0:
        N = int(l)

    il += 1

# Output PDB
fout = open(args.out, 'w')
fout.write('%i\n' % N)
fout.write('\n')

# Open DCD and read the header
dcd = DcdFile(args.dcd)
dcd.open_to_read()
dcd.read_header()

num_dcd_frames = dcd.count_frame()

if args.frame is None:
    args.frame = num_dcd_frames - 1

# skip
skipped = dcd.skip_as_many_as_possible_upto(args.frame)
if skipped < args.frame:
    print('Only %i frames could be read from the dcd file, which is less than the one specified as --frame %i.' % (skipped, args.frame_begin))
    print('No output file generated.')
    sys.exit(2)


# read and write
if not dcd.has_more_data() :
    print('The number of frames is invalid.')
    print('Header information:')
    dcd.show_header()
    dcd.close()
    fout.close()
    sys.exit(2)
    
struct = dcd.read_onestep()

""" Move to the origin """
if args.flg_origin or args.flg_4500:
    com = [0.0, 0.0, 0.0]
    for v in struct:
        com = [com[i]+v[i] for i in range(3)]

    com = [com[i]/float(len(struct)) for i in range(3)]

    for i in range(len(struct)):
        if args.flg_4500:
            struct[i] = [struct[i][j] - com[j] + 4500.0 for j in range(3)]
        else:
            struct[i] = [struct[i][j] - com[j] for j in range(3)]

for i in range(N):
    fout.write(seq[i])
    fout.write("  %5.2f  %5.2f  %5.2f\n" % (struct[i][0], struct[i][1], struct[i][2]))

fout.close()
