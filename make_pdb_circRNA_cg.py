#!/usr/bin/env python

import math
import argparse

parser = argparse.ArgumentParser(
             formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('outfile', type=argparse.FileType('w'), help='output filename')

group_seq = parser.add_mutually_exclusive_group(required=True)
group_seq.add_argument('--seq', help='Sequence')
group_seq.add_argument('--seqfile', type=argparse.FileType('r'), help='Sequence FASTA file')

parser.add_argument('--b', type=float, default=5.5, help='Distance between phosphates of adjacent nucleotides')
parser.add_argument('--l', type=float, default=4.2, help='Distance between Phosphate and Sugar')

args = parser.parse_args()

if args.seq is not None:
    seq = args.seq
    n_nt = len(seq)

elif args.seqfile is not None:
    seq = ''
    for l in args.seqfile:
        if l[0] == '>' or l[0] == '#' or l[0] == ';':
            continue
        seq += l.strip()
    n_nt = len(seq)

n_nt = len(seq)
print ('Sequence: %s' % (seq,))
print ('#nucleotides: %i' % (n_nt,))
b = args.b   #5.5  # Distance between P and P
l = args.l   #4.2  # Distance between P and S

theta = 2 * math.pi / float(n_nt)
r = b / theta

for i in range(n_nt):
    s = seq[i-1]

    x = r * math.cos(theta * i)
    y = 0
    z = r * math.sin(theta * i)
    args.outfile.write('ATOM  %5i  P   R%s  A %3i    %8.3f%8.3f%8.3f\n' % (3*i+1, s, i+1, x,y,z))

    x = r * math.cos(theta * (i+0.5))
    y = - math.sqrt(l*l - b*b*0.25) 
    z = r * math.sin(theta * (i+0.5))
    args.outfile.write('ATOM  %5i  S   R%s  A %3i    %8.3f%8.3f%8.3f\n' % (3*i+2, s, i+1, x,y,z))

    x = r * math.cos(theta * (i+0.5))
    y = - math.sqrt(l*l - b*b*0.25) - l
    z = r * math.sin(theta * (i+0.5))
    args.outfile.write('ATOM  %5i  %sb  R%s  A %3i    %8.3f%8.3f%8.3f\n' % (3*i+3, s, s, i+1, x,y,z))

args.outfile.close()
