#!/usr/bin/env python
'''
@author: Naoto Hori
'''

import sys
import argparse
from lop.file_io.dcd import DcdFile

parser = argparse.ArgumentParser(description='Script to check periodic boundary box information in DCD file',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--frame', default=0,
                    action='store', type=int, help='Frame number (0 start)')
parser.add_argument('dcdfile', help='Input DCD file')

args = parser.parse_args()

if args.frame < 0:
    print('Error: --frame is invalid.')
    sys.exit(2)

dcd = DcdFile(args.dcdfile)
dcd.open_to_read()
dcd.read_header()

for i in range(args.frame):
    dcd.skip_onestep()

    if not dcd.has_more_data() :
        print('Error: --frame is invalid.')
        print('Header information:')
        dcd.show_header()
        dcd.close()
        sys.exit(2)

dcd.read_onestep()

dcd.show_unit_cell()
dcd.close()
