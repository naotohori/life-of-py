#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2013/08/05
@author: Naoto Hori
'''

import sys

if len(sys.argv) != 3:
    print 'Usage SCRIPT [input pdb] [output prefix]'
    sys.exit(2)
    
f_in = open(sys.argv[1])
pfx = sys.argv[2]

flg_open = False
for l in f_in:
    if l[0:5] == 'MODEL':
        filename = pfx + ('%08i' % int(l[5:])) + '.pdb'
        f_out = open(filename,'w')
        flg_open = True
    elif (l[0:6] == 'REMARK' or l[0:4] == 'ATOM'   or
          l[0:6] == 'HETATM' or l[0:3] == 'TER'    ):
        if flg_open:
            f_out.write(l)
        else:
            print('file is not open for a line:')
            print('     '+l)
    elif l[0:6] == 'ENDMDL':
        f_out.write('END')
        f_out.close()
        flg_open = False
        
        
        
