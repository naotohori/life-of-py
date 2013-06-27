#!/usr/bin/env python
'''
Created on 2011/12/01
@author: Naoto Hori
'''

import sys

if len(sys.argv) != 3:
    print ('Usage: SCRIPT [IN opt] [OUT opt.ene]')
    sys.exit(2)

f_out = open(sys.argv[2], 'w')

first = -1
s = 0.0
for line in open(sys.argv[1],'r'):
    if line.find('#') != -1:
        continue
    linesp = line.split()
    i = int(linesp[0])
    if first < 0:
        first = i
    elif i == first:
        f_out.write(' %15.6f\n' % (s,))
        s = 0.0
    else:
        s += float(linesp[-1])






