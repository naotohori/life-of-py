#!/usr/bin/env python
# 2011/03/17 coded by Naoto HORI

import sys
from cafysis.file_io.dcd import DcdFile
from numpy import zeros, float32

if len(sys.argv) != 3:
    print('Usage: % SCRIPT [input DCD] [output]')
    sys.exit(2)
    
filename_out = sys.argv[2]

dcd = DcdFile(sys.argv[1])
dcd.open_to_read()
dcd.read_header()
nmp = dcd.get_header().nmp_real

print('nmp=%i\n'%nmp)
num_model = 0
covariance = zeros((nmp*3, nmp*3), dtype=float32)
variance = zeros((nmp*3,), dtype=float32)

while dcd.has_more_data() :
    data = dcd.read_onestep()
    num_model += 1
    print(num_model)
    
#    for i_mp in xrange(nmp) :
#        for i_xyz in xrange(3) :
#            i = i_mp + i_xyz
#            for j_mp in xrange(i_mp+1) :
#                for j_xyz in xrange(3):
#                    j = j_mp + j_xyz
#                    covariance[i,j] += data[i_mp][i_xyz] * data[j_mp][j_xyz]
#            variance[i] += data[i_mp][i_xyz]
#            
#a = zeros((nmp*3,nmp*3), dtype=float32)
#for i in xrange(nmp):
#    for j in xrange(i+1) :
#        a[i,j] = covariance[i,j] / float(num_model) \
#                - variance[i] * variance[j] / (float(num_model) **2)
#    for j in xrange(i+1, nmp) :
#        a[i,j] = covariance[j,i] / float(num_model) \
#                - variance[i] * variance[j] / (float(num_model) **2)

file_out = open(filename_out,'w')
file_out.write('%i\n' % nmp*3)
for i in range(nmp*3) :
    for j in range(i,nmp*3):
        file_out.write('%22.20f\n' % a[i,j])

file_out.close()
