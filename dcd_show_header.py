#!/usr/bin/env python
'''
@author: Naoto Hori
'''

import sys
from lop.file_io.dcd import DcdFile

if (len(sys.argv) != 2):
    print('Usage: % SCRIPT [filename]')
    sys.exit(2)
    
dcd = DcdFile(sys.argv[1])
dcd.open_to_read()
dcd.read_header()
dcd.show_header()
dcd.close()
