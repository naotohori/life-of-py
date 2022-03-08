#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@author: Naoto Hori
'''

import sys
from lop.file_io.pdb import PdbFile

if len(sys.argv) < 5:
    print('Usage: % SCRIPT [average file] [input pdb] [ev file] [,ev file ....] [output PC file]')
    sys.exit(2)
    
data_ave = []
for line in open(sys.argv[1], 'r'):
    data_ave.append(float(line.strip()))

print(("#Dimension: %s" % len(data_ave)))

f_pdb = PdbFile(sys.argv[2])
f_pdb.open_to_read()

f_out = open(sys.argv[-1], 'w')

# Read eigen values
num_ev = len(sys.argv) - 4
ev = []
for i in range(num_ev) :
    f_ev = open(sys.argv[i+3], 'r')
    ev_tmp = [] 
    for line in f_ev :
        if line.find('#') != -1: continue
        ev_tmp.append(float(line.strip()))
    ev.append(ev_tmp)
    f_ev.close()
    
# Check for ev
num_dimension = len(ev[0])
for v in ev :
    if len(v) != num_dimension :
        print('len(v) != num_dimension, %i' % num_dimension)
        sys.exit(2)

#debug
#for i in xrange(num_ev):
#    print '# ev %i' % i
#    for value in ev[i] :
#        print value

# Calculate principal coordinate
chains = f_pdb.read_all()

for v in ev:
    x = 0.0
    idx = 0
    for c in chains:
        for i in range(c.num_atom()):
            c.get_atom(i).xyz.x
            x += (c.get_atom(i).xyz.x - data_ave[idx]) * v[idx]
            idx += 1
            x += (c.get_atom(i).xyz.y - data_ave[idx]) * v[idx]
            idx += 1
            x += (c.get_atom(i).xyz.z - data_ave[idx]) * v[idx]
            idx += 1
    f_out.write('%12.5f ' % x)
f_out.write('\n')
            
f_pdb.close()
f_out.close()

