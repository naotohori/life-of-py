#!/usr/bin/env python

import sys

f_out = open('new.fasta', 'w')

def output(seq, f):
    s = len(seq)
    while s > 0:
        for i in range(5):
            if len(seq) > 10:
                f.write(seq[0:10])
                seq = seq[10:]
                s -= 10
                f.write(' ')
            else:
                f.write(seq)
                s = 0
                break
        f.write('\n')
    f.write('\n')

seq = ''
flg = False
for l in open(sys.argv[1]):

    if l.startswith('>'):
        if flg:
            output(seq, f_out)
            flg = False
        f_out.write(l)

    elif len(l.strip()) == 0:
        continue

    else:
        if flg:
            seq += l.strip()
        else:
            flg = True
            seq = l.strip()

if flg:
    output(seq, f_out)
    flg = False

f_out.close()
