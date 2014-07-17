#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2014/06/17
@author: Naoto Hori
'''

import sys
from cafysis.file_io.dcd import DcdFile

if len(sys.argv) != 2:
    print ('Usage: SCRIPT [dcd file]')
    sys.exit(2)

dcd = DcdFile(sys.argv[1])
dcd.open_to_read()

dcd.read_header()

def error_no_data() :
    print 'The number of frames is invalid.'
    print 'Header information:'
    dcd.show_header()
    sys.exit(2)


icount = 0
while dcd.has_more_data():
    dcd.skip(1)
    icount += 1

dcd.close()

print ('# frames = ', icount)
