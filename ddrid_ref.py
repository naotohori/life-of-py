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
    print 'Usage: % SCRIPT [input drid] [input drid1] [output dDRID]'
    sys.exit(2)

drid = DridFile(sys.argv[1])
drid.open_to_read()
drid.read_header()

drid1 = DridFile(sys.argv[2])
drid1.open_to_read()
drid1.read_header()

#if len(drid.get_header().centroids) != len(drid1.get_header().centroids):
if not np.array_equal(drid.get_header().mask, drid1.get_header().mask):
    print 'Two masks are not identical'
    sys.exit(2)

x_ref = drid1.read_onestep()
drid1.close()

f_out = open(sys.argv[-1], 'w')
while drid.has_more_data():
    x = drid.read_onestep()
    ddrid = py_ddrid.ddrid(x, x_ref)

    f_out.write('%f\n' % (ddrid,))

drid.close()
f_out.close()
