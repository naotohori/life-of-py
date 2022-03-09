#!/usr/bin/env python

from lop.file_io.sisbp import SisbpFile
import sys

f = SisbpFile(sys.argv[1])
f.open_to_read()

while f.has_more_data():

    pairs, energies = f.read_onestep()

    for p,e in zip(pairs, energies):
        print(f" ({p[0]:d}, {p[1]:d}: {e:f})", end='')
    print("")

f.close()
