#!/usr/bin/env python

import sys
import scipy.cluster
import numpy as np
from CalcRMSD import calcrmsd
from lop.file_io.dcd import DcdFile

if len(sys.argv) != 4:
    print('Usage: % SCRIPT [input DCD] [#frame to skip] [output prefix]')
    sys.exit(2)

nskip = int(sys.argv[2])
prefix = sys.argv[-1]

dcd = DcdFile(sys.argv[1])
dcd.open_to_read()
dcd.read_header()

dist_array = []

while dcd.has_more_data():
    ref = dcd.read_onestep_np()
    dcd.set_mark()

    if dcd.has_more_data():
        while dcd.has_more_data():
            data = dcd.read_onestep_np()
            dist_array.append( calcrmsd(ref.T, data.T) )
        dcd.go_mark()
    else:
        break

dcd.close()

print(len(dist_array))

#z = scipy.cluster.hierarchy.ward(dist)  # This does not work for some reason 
z = scipy.cluster.hierarchy.linkage(dist_array, method='ward', metric='euclidean')
## Other methods
#z = scipy.cluster.hierarchy.linkage(data, method='single', metric='euclidean')
#z = scipy.cluster.hierarchy.linkage(data, method='average', metric='euclidean')

f_out = open(prefix+'.dcd.cls.z', 'w')
for iz, zz in enumerate(z):
    f_out.write("%5i %5i %5i %f %i\n" % (iz+1, zz[0],zz[1],zz[2],zz[3]))
f_out.close()


f_out = open(prefix+'.dcd.cls.ncls', 'w')
dmin = z[0][2]
dmax = z[-1][2]
for i in range(100):
    d = dmin + i * (dmax - dmin) / 100.0
    ncls = scipy.cluster.hierarchy.fcluster(z, d, criterion='distance').max()
    f_out.write('%f %i\n' % (d, ncls))
f_out.close()
