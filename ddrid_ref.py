#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Created on 2015/05/26
@author: Naoto Hori
'''

import sys
import numpy as np
from cafysis.file_io.drid import DridFile, DridHeader
from cafysis.lib_f2py import py_ddrid

if len(sys.argv) != 4:
    print('Usage: % SCRIPT [input drid] [reference drid] [output dDRID]')
    sys.exit(2)

drid = DridFile(sys.argv[1])
drid.open_to_read()
drid.read_header()

drid_ref = DridFile(sys.argv[2])
drid_ref.open_to_read()
drid_ref.read_header()

if not np.array_equal(drid.get_header().mask, drid_ref.get_header().mask):
    print('Two masks are not identical')
    sys.exit(2)

x_ref = drid_ref.read_onestep()
drid_ref.close()

f_out = open(sys.argv[-1], 'w')
while drid.has_more_data():
    x = drid.read_onestep()
    ddrid = py_ddrid.ddrid(x, x_ref)

    f_out.write('%f\n' % (ddrid,))

drid.close()
f_out.close()
