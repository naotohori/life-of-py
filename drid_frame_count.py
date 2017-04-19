#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2017/04/19
@author: Naoto Hori
'''

import sys
from cafysis.file_io.drid import DridFile

def count(path):

    drid = DridFile(path)
    drid.open_to_read()

    drid.read_header()

    icount = 0
    while drid.has_more_data():
        drid.skip(1)
        icount += 1

    drid.close()

    return icount


if __name__ == "__main__":

    if len(sys.argv) != 2:
        print ('Usage: SCRIPT [drid file]')
        sys.exit(2)

    print ('# frames = ', count(sys.argv[1]))
