#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@author: Naoto Hori
'''

import sys
import argparse

from lop.file_io.pdb import PdbFile
from lop.rna_ss_convert import convert_format, write_output_file
from lop.file_io.sisbp import SisbpFile

parser = argparse.ArgumentParser(description='Extract sisbp frame as various format', 
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--frame', dest='frame', default=None,
                    action='store', type=int, help='Frame ID (starting from 0, default: last frame)')

group_seq = parser.add_mutually_exclusive_group(required=True)
group_seq.add_argument('--pdb', type=PdbFile, help='PDB file')
group_seq.add_argument('--seqtext', help='Sequence')
group_seq.add_argument('--seqfile', type=argparse.FileType('r'), help='Sequence FASTA file')

parser.add_argument('bpfile', type=SisbpFile, help='input sisbp file')
parser.add_argument('outfile', help='Output file with appropriate extension (.db/.ct/.bpseq)')

args = parser.parse_args()

if args.frame is not None:
    if args.frame < 0:
        print('Error: Frame ID must be >= 0')
        sys.exit(2)

if args.pdb is not None:
    #pdb = PdbFile(sys.argv[1])
    args.pdb.open_to_read()
    chains = args.pdb.read_all()
    args.pdb.close()

    seq = ['',]
    n_nt = [0,]

    for c in chains:
        n_nt.append(c.num_res())

        s = []
        for r in c.residues:
            # "RA " ---> "A"
            s.append(r.atoms[0].res_name.strip()[1])

        seq.append(s)

elif args.seqtext is not None:
    # This case, multiple chains are not allowed.
    seq = ['', args.seqtext,]
    n_nt = [0, len(seq[1]),]

elif args.seqfile is not None:
    seq = ['',]
    n_nt = [0,]
    s = ''
    flg_reading = False

    for l in args.seqfile:
        if l[0] == '>' or l[0] == '#' or l[0] == ';':
            if flg_reading:
                seq.append(s)
                n_nt.append(len(s))
                flg_reading = False
                s = ''
            continue

        elif len(l.strip()) == 0:
            continue

        else:
            flg_reading = True
            s += l.strip()

    if flg_reading:
        seq.append(s)
        n_nt.append(len(s))

else:
    print ('Error: this error should be detected by mutually_exclusive_group of the parser')
    sys.exit(2)

f = args.bpfile
f.open_to_read()

if args.frame is None:
    frame_skip = f.count_frame() - 1
else:
    frame_skip = args.frame

fout = open(args.outfile, 'w')

# skip
skipped = f.skip_as_many_as_possible_upto(frame_skip)
if skipped < frame_skip:
    print('Only %i frames could be read from the sisbp file, which is less than the one specified as --frame %i.' % (skipped, args.frame_begin))
    print('No output file generated.')
    sys.exit(2)

# read and write
if not f.has_more_data() :
    print('Only %i frames could be read from the sisbp file, which is less than the one specified as --frame %i.' % (skipped + 1, args.frame_begin))
    print('No output file generated.')
    f.close()
    sys.exit(2)

pairs, energies = f.read_onestep()

f.close()

write_output_file(seq[1], pairs, args.outfile)
