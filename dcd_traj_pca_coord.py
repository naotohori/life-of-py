#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2011/07/30
@author: Naoto Hori
'''

import sys
from lop.file_io.dcd import DcdFile
from lop.file_io.pdb import PdbFile


if len(sys.argv) != 5:
    print(' Usage: % SCRIPT [input DCD] [average PDB] [input eigen vector] [output] ')
    sys.exit(2)
    
dcd = DcdFile(sys.argv[1])
dcd.open_to_read()
dcd.read_header()
nmp = dcd.get_header().nmp_real

# reading average
f_pdb = PdbFile(sys.argv[2])
f_pdb.open_to_read()
chains = f_pdb.read_all()
f_pdb.close()
ave = []
for c in chains:
    imp = 0
    for i in range(c.num_atom()) :
        imp += 1
        ave.append(c.get_atom(i).xyz.x)
        ave.append(c.get_atom(i).xyz.y)
        ave.append(c.get_atom(i).xyz.z)
if len(ave) != 3*nmp :
    print('Error: len(ave) != 3*nmp')
    sys.exit(2)

# reading pca eigenvector
f_pca = open(sys.argv[3], 'r')
ev = []
for line in f_pca:
    if line.find('#') != -1:
        continue
    ev.append(float(line))
f_pca.close()
print(ev)
if len(ev) != 3*nmp :
    print('Error: len(ev) != 3*nmp')
    sys.exit(2)

# loop for DCD
f_out = open(sys.argv[4], 'w')
nframe = 0
while dcd.has_more_data() :
    data = dcd.read_onestep()
    nframe += 1
    if nframe % 1000 == 0:
        print(nframe)
    coord = 0.0
    for i in range(nmp) :
        coord += ( (data[i][0] - ave[3*i+0]) * ev[3*i+0] 
                 + (data[i][1] - ave[3*i+1]) * ev[3*i+1] 
                 + (data[i][2] - ave[3*i+2]) * ev[3*i+2] )
    f_out.write('%12.5f\n' % coord)
dcd.close()
f_out.close()

##!/usr/bin/env python
## 2011/07/30 coded by Naoto HORI
#
#import sys
#from file_dcd import DcdFile
#
#if len(sys.argv) != 4:
#    print ' Usage: % SCRIPT [input DCD] [input eigen vector] [output] '
#    sys.exit(2)
#    
#dcd = DcdFile(sys.argv[1])
#dcd.open_to_read()
#dcd.read_header()
#nmp = dcd.get_header().nmp_real
#
#f_pca = open(sys.argv[2], 'r')
#ev = []
#for line in f_pca:
#    if line.find('#') != -1:
#        continue
#    ev.append(float(line))
#f_pca.close()
#print ev
#if len(ev) != 3*nmp :
#    print 'Error: len(ev) != 3*nmp'
#    sys.exit(2)
#
#f_out = open(sys.argv[3], 'w')
#nframe = 0
#while dcd.has_more_data() :
#    data = dcd.read_onestep()
#    nframe += 1
#    if nframe % 1000 == 0:
#        print nframe
#    coord = 0.0
#    for i in xrange(nmp) :
#        coord += ( data[i][0] * ev[3*i+0] 
#                 + data[i][1] * ev[3*i+1] 
#                 + data[i][2] * ev[3*i+2] )
#    f_out.write('%12.5f\n' % coord)
#dcd.close()
#f_out.close()