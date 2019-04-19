#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2011/05/05
@author: Naoto Hori
'''

'''
CafeMol出力のPDBファイルを、CafeMolのinputファイルの"pdb"行を元に
chainごとに分割する。
'''

import sys
import re

if len(sys.argv) != 4:
    print()
    print('Usage: % SCRIPT [cafemol pdb file] [cafemol input file]  [output prefix]')
    print()
    sys.exit(2)

f_pdb = open(sys.argv[1], 'r')
f_inp = open(sys.argv[2], 'r')
out_prefix = sys.argv[3]

flg_to_read = False
re_filename = re.compile('(\d*)\s*(protein|rna)\s*(\S*)\s*') 
data = {}
for line in f_inp :
    if re.match('<<<< unit_and_state', line) :
        flg_to_read = True
        continue
    if flg_to_read :
        if line[0:1] == '*' :
            continue
        if re.match('>>>>', line) :
            break
        mo = re_filename.match(line)
        if mo :
            i = int( mo.groups()[0] )
            data[i] = (mo.groups()[1:])
f_inp.close()

for i in range(1, len(data)+1) :
    f_out = open(out_prefix+data[i][1], 'w')
    for line in f_pdb :
        if line[:2] == '>>' :
            break
        if line[:4] == 'ATOM' :
            f_out.write(line)
    f_out.close()
