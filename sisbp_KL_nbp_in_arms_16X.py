#!/usr/bin/env python

from lop.file_io.sisbp import SisbpFile
import sys
import argparse

NMOL = 16
NMP_PER_MOL = 63 * 4
NARM = 4

MPS_IN_ARMS = [[10, 11, 12, 13, 14, 15], [73, 74, 75, 76, 77, 78], [136, 137, 138, 139, 140, 141], [199, 200, 201, 202, 203, 204]]

parser = argparse.ArgumentParser(
         description='Calculate number of base pairs from bp file',
         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('bpfile', type=SisbpFile, help='input bp file')
parser.add_argument('outfile', default='nbp.out', help='output filename')

args = parser.parse_args()

f = args.bpfile

f.open_to_read()
fout = open(args.outfile, 'w')

while f.has_more_data():

    pairs, energies = f.read_onestep()

    mp_forming_bp = set()
    for imp, jmp in pairs:
        mp_forming_bp.add(imp)
        mp_forming_bp.add(jmp)

    for imol in range(NMOL):
        for iarm in range(NARM):
            n = 0
            for imp in MPS_IN_ARMS[iarm]:
                imp_local = imp + 1 + NMP_PER_MOL * imol
                if imp_local in mp_forming_bp:
                    n += 1
            fout.write(f' {n}')
    fout.write('\n')

f.close()
fout.close()
