#!/usr/bin/env python

import sys

if len(sys.argv) != 3:
    print('Usage: SCRIPT [DCD] [output matrix file]')
    sys.exit(2)

from cafysis.file_io.dcd import DcdFile
from CalcRMSD import calcrmsd

dcd = DcdFile(sys.argv[1])
dcd.open_to_read()
dcd.read_header()

f_out = open(sys.argv[-1],'w')

i = -1
dcd.set_mark()
while dcd.has_more_data():
    i += 1
    ref = dcd.read_onestep_np()
    dcd.set_mark()

    j = i
    while dcd.has_more_data():
        j += 1
        data = dcd.read_onestep_np()
        rmsd = calcrmsd(ref.T, data.T)
        f_out.write('%i %i %f\n' % (i,j,rmsd))

    if j==i:
        break
    else:
        dcd.go_mark()

dcd.close()
f_out.close()
