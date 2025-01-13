#!/usr/bin/env python
'''
@author: Naoto Hori
'''

import sys
import argparse
from lop.file_io.dcd import DcdFile

parser = argparse.ArgumentParser(description='Script to check periodic boundary box information in DCD file',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
group_frame = parser.add_mutually_exclusive_group(required=True)
group_frame.add_argument('--frame', action='store', type=int, help='Frame number (0 start)')
group_frame.add_argument('--allframes', default=False,
                    action='store_true', help='Show for all frames.')
parser.add_argument('dcdfile', help='Input DCD file')

args = parser.parse_args()

dcd = DcdFile(args.dcdfile)
dcd.open_to_read()
dcd.read_header()

if args.allframes:
    i = 0
    while dcd.has_more_data():

        dcd.read_onestep()
        print (f'{i} ', end='')
        dcd.show_unit_cell_oneline()
        i += 1

else:

    if args.frame < 0:
        print('Error: --frame is invalid.')
        sys.exit(2)

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
