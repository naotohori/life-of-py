#!/usr/bin/env python

from cafysis.file_io.drid import DridFile
import numpy as np
import sys
import math

if len(sys.argv) != 3:
    print('Usage: % SCRIPT [input DRID] [output pairwise]')
    sys.exit(2)

f_out = open(sys.argv[-1],'w')

drid = DridFile(sys.argv[1])
drid.open_to_read()
drid.read_header()

ndim = 3*drid._header.n_centroids()
nframe = 0
while drid.has_more_data():
    drid.read_onestep()
    nframe += 1

drid.rewind()

data = np.zeros( (nframe, ndim) )
k = 0
while drid.has_more_data():
    data[k,:] = drid.read_onestep()
    k += 1

drid.close()

for i in range(nframe):
    for j in range(i+1, nframe):
        f_out.write('%i %i %f\n' % (i,j, math.sqrt( np.square(data[i,:] - data[j,:]).sum() / float(ndim)) ))

f_out.close()


