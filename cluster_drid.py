#!/usr/bin/env python

import sys
import scipy
import scipy.spatial
import scipy.cluster
import numpy as np
from lop.file_io.drid import DridFile
import math

if len(sys.argv) != 4:
    print('Usage: % SCRIPT [input DRID] [#frame to skip] [output prefix]')
    sys.exit(2)

nskip = int(sys.argv[2])
prefix = sys.argv[-1]

drid = DridFile(sys.argv[1])
drid.open_to_read()
drid.read_header()

nc = 3*drid.get_header().n_centroids()
nframe = 0
while drid.has_more_data():
    drid.read_onestep()
    nframe += 1
    if drid.has_more_data():
        drid.skip(nskip-1)

drid.rewind()

data = np.zeros( (nframe, nc) )
k = 0
while drid.has_more_data():
    data[k,:] = drid.read_onestep()
    if drid.has_more_data():
        drid.skip(nskip-1)
    k += 1

drid.close()


###! DRID has to be compared by RMSD(DRID)
###! Becauase the linkage method computes only "euclidean" distance.
###! division factor sqrt(nc)  (nc = 3 * n_centroids) should be included in the data here.
data = data / math.sqrt(nc)

#z = scipy.cluster.hierarchy.ward(dist)  # This does not work for some reason 
z = scipy.cluster.hierarchy.linkage(data, method='ward', metric='euclidean')
## Other methods
#z = scipy.cluster.hierarchy.linkage(data, method='single', metric='euclidean')
#z = scipy.cluster.hierarchy.linkage(data, method='average', metric='euclidean')

f_out = open(prefix+'.drid.cls.z', 'w')
for iz, zz in enumerate(z):
    f_out.write("%5i %5i %5i %f %i\n" % (iz+1, zz[0],zz[1],zz[2],zz[3]))
f_out.close()


f_out = open(prefix+'.drid.cls.ncls', 'w')
dmin = z[0][2]
dmax = z[-1][2]
for i in range(100):
    d = dmin + i * (dmax - dmin) / 100.0
    ncls = scipy.cluster.hierarchy.fcluster(z, d, criterion='distance').max()
    f_out.write('%f %i\n' % (d, ncls))
f_out.close()
