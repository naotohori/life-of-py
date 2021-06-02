#!/usr/bin/env python

from cafysis.file_io.sisbp import SisbpFile
import sys

f = SisbpFile(sys.argv[1])
f.open_to_read()

while f.has_more_data():

    pairs = f.read_onestep()

    for p in pairs:
        print(f" ({p[0]:d}, {p[1]:d})", end='')
    print("")

f.close()
