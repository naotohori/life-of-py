#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2011/03/17 (dcd_frame_extract.py)
Added stride mode on 2011/07/29 (dcd_frame_extract.py)
Modified from dcd_frame_extract.py on 2011/12/20
@author: Naoto Hori
'''

import sys
import argparse

from lop.file_io.dcd import DcdFile
from lop.file_io.pdb import PdbFile

parser = argparse.ArgumentParser(description='Convert DCD to Movie (sequence of pdb format)',
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

#parser.add_argument('--fit', dest='flg_fit', default=False,
#                    action='store_true', help='Fit to the reference PDB')

#parser.add_argument('--id-digit', dest='id_digit', default=0,
#                    action='store', type=int, help='Number of digit of ID for MODEL')

parser.add_argument('dcd', help='Input DCD file')
parser.add_argument('pdb', help='Reference PDB file')
parser.add_argument('out', help='Output movie file')

args = parser.parse_args()

#if (not len(sys.argv) in (6, 7)):
#    print('Usage: % SCRIPT [input DCD] [beginning (0)] [end] [reference PDB] [output movie]')
#    print('       % SCRIPT [input DCD] [beginning (0)] [end] [stride] [reference PDB] [output movie]')
#    sys.exit(2)

#if args.flg_fit:
#    from Superimpose import superimpose

if args.frame_begin < 0:
    print('Error: Beginning frame must be > 0')
    sys.exit(2)

if args.frame_stride <= 0:
    print('Error: Frame stride must be > 0')
    sys.exit(2)

# Read the reference PDB
pdb = PdbFile(args.pdb)
pdb.open_to_read()
chains = pdb.read_all()
pdb.close()

# Output PDB
movie = PdbFile(args.out)
movie.open_to_write()

# Open DCD and read the header
dcd = DcdFile(args.dcd)
dcd.open_to_read()
dcd.read_header()

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
        movie.close()
        sys.exit(2)
        
    if iframe % args.frame_stride != 0 :
        dcd.skip(1)
        i_orig += 1
        continue
    
    struct = dcd.read_onestep()

    """ Move to the origin """
    if args.flg_origin:
        com = [0.0, 0.0, 0.0]
        for v in struct:
            com = [com[i]+v[i] for i in range(3)]

        com = [com[i]/float(len(struct)) for i in range(3)]

        for i in range(len(struct)):
            struct[i] = [struct[i][j] - com[j] for j in range(3)]

    iatom = 0
    for c in chains:
        for r in c.residues:
            for a in r.atoms:
                a.xyz.put_as_list(struct[iatom])
                iatom += 1

    movie.modelID = iframe + args.id_offset
    movie.set_remark('ORIGINAL_FILE %s' % (args.dcd,))
    movie.set_remark('ORIGINAL_FRAME %i' % (i_orig,))
    movie.write_all(chains)
    i_orig += 1

movie.close()
pdb.close()
