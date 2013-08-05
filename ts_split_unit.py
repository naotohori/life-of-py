#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2013/08/01
@author: Naoto Hori
'''

import sys

if len(sys.argv) != 2:
    print ' Usage: % SCRIPT [ts file]'
    sys.exit(2)
    
tsfilename = sys.argv[1]
f_ts = open(tsfilename,'r')

## Count the number of units
# Read to the first 'all'
for l in f_ts:
    if l.find('#all') != -1:
        break
# Read to the next 'all'
u_list = []
u_list.append('all')
for l in f_ts:
    if l.find('#all') != -1:
        break
    if l.find('#') != -1:
        u_list.append(l.split()[0][1:])

## Generate output files
f_list = []
for u in u_list:
    f_list.append(open(tsfilename+'.'+u,'w'))
    
## Main
f_ts.seek(0)
for l in f_ts:
    for u,f in zip(u_list,f_list):
        if l.find('#'+u) == 0:
            f.write(l[len(u)+1:])
            
f_ts.close()
for f in f_list:
    f.close()
    