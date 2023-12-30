#!/usr/bin/env python

from lop.file_io.sisbp import SisbpFile
import sys

if len(sys.argv) != 2:
    print('Usage: SCRIPT [sisbp file]')
    sys.exit(2)

f = SisbpFile(sys.argv[1])
f.open_to_read()
print(f.count_frame())
f.close()
