#!/usr/bin/env python

import sys
import math
import numpy as np
from Superimpose import superimpose
import scipy.cluster.hierarchy
from lop.file_io.dcd import DcdFile

if len(sys.argv) != 6:
    print('Usage: SCRIPT [DCD file] [Native (reference) DCD] [prefix] [cutoff] [nskip (to calculate frame id)]')
    sys.exit(2)

dcd_filepath = sys.argv[1]
dcd_ref_filepath = sys.argv[2]
prefix = sys.argv[3]
cutoff_char = sys.argv[4]
cutoff = float(sys.argv[4])
nskip = int(sys.argv[5])

z = []
for l in open(prefix+'.cls.z'):
    lsp = l.split()
    z.append( [int(lsp[1]), int(lsp[2]), float(lsp[3]), int(lsp[4])] )


fcls = scipy.cluster.hierarchy.fcluster(z, cutoff, criterion='distance')   # 300

#print len(fcls)
#print fcls
f_out = open(prefix+'.cls_%s' % cutoff_char, 'w')
for f in fcls:
    f_out.write('%i\n' % (f,))
f_out.close()
#
##L, M = scipy.cluster.hierarchy.leaders(z, fcls)
##print 'L:',len(L)
##print L
##print 'M:',len(M)
##print M
#

ncls = max(fcls)


f_out = open(prefix+".cls_%s.summary" % cutoff_char, 'w')
f_out.write('#')
for arg in sys.argv:
    f_out.write(' %s' % (arg,))
f_out.write('\n')

# Size of cluster
f_out.write('#Size of cluster\n')
for i in range(1,ncls+1):
    f_out.write('%i %i\n' % (i, fcls.tolist().count(i)))
f_out.write('\n\n')

# Member of cluster
f_out.write('### Members of cluster\n')
for i in range(1,ncls+1):
    f_out.write('#cluster %i\n' % (i,))
    for iframe, f in enumerate(fcls): 
        if f == i:
            f_out.write("%i %i\n" % (iframe, iframe*nskip))
    f_out.write('\n\n')
f_out.close()



############################# Center
dcd_ref = DcdFile(dcd_ref_filepath)
dcd_ref.open_to_read()
dcd_ref.read_header()
data_ref = dcd_ref.read_onestep_np()
dcd_ref.close()

dcd = DcdFile(dcd_filepath)
dcd.open_to_read()
dcd.read_header()

nmp = dcd._header.nmp_real
cls_centroids = []
for icls in range(ncls):
    cls_centroids.append( np.zeros((nmp,3)) )

k = 0
while dcd.has_more_data():
    data = dcd.read_onestep_np()
    rmsd = superimpose(data_ref.T, data.T)

    cls_centroids[fcls[k]-1] += data[:,:]
    k+=1

for icls in range(ncls):
    #cls_centroids[icls] /= float(k)
    cls_centroids[icls] /= float( fcls.tolist().count(icls+1) )

#cls_nearest_DCD = []
cls_nearest_node = []
cls_nearest_RMSD = []
cls_average_RMSD = []
cls_num_node = []
for icls in range(ncls):
    #cls_nearest_DCD.append( np.zeros((nmp,3)) )
    cls_nearest_node.append(-1)
    cls_nearest_RMSD.append(9999999.)
    cls_average_RMSD.append(0.0)
    cls_num_node.append(0)

dcd.rewind()
k = 0
while dcd.has_more_data():
    data = dcd.read_onestep_np()

    icls = fcls[k] - 1
    rmsd = superimpose(cls_centroids[icls].T, data.T)

    if rmsd < cls_nearest_RMSD[icls]:
        cls_nearest_RMSD[icls] = rmsd
        #cls_nearest_DCD[icls] = data
        #cls_nearest_node[icls] = k + 1
        cls_nearest_node[icls] = k
    cls_average_RMSD[icls] += rmsd
    cls_num_node[icls] += 1
    k += 1

# Check this code correctly working
for icls in range(ncls):
    if cls_num_node[icls] != fcls.tolist().count(icls+1):
        print('Error: cls_num_node[icls] != fcls.count(icls+1)')
        sys.exit(2)

for icls in range(ncls):
    cls_average_RMSD[icls] /= float(cls_num_node[icls])

f_out = open('%s.cls_%s.centroids' % (prefix, cutoff_char),'w')
f_out.write('#clsID node  RMSD  RMSD/sqrt(nmp) <RMSD> <RMSD>/sqrt(nmp)\n')
for icls in range(ncls):
    f_out.write('%i %i %5.2f %f %5.2f %f\n' % (icls+1, cls_nearest_node[icls], 
            cls_nearest_RMSD[icls], cls_nearest_RMSD[icls]/math.sqrt(nmp),
            cls_average_RMSD[icls], cls_average_RMSD[icls]/math.sqrt(nmp)) )
f_out.write('\n')


dcd_out = DcdFile('%s.cls_%s.centroids.dcd' % (prefix, cutoff_char))
dcd_out.open_to_write()
dcd_out.set_header( dcd.get_header() )
dcd_out.write_header()

for icls in range(ncls):
    dcd_out.write_onestep( cls_centroids[icls] )
dcd_out.close()

#f_out.write('#centroid DCD\n')
#for icls in range(ncls):
#    f_out.write('%i' % (icls+1,))
#    for c in cls_centroids[icls]:
#        f_out.write(' %f' % (c,))
#    f_out.write('\n')
#f_out.write('\n')
#
#f_out.write('#nearest DCD\n')
##for icls in range(ncls):
#    f_out.write('%i' % (icls+1,))
#    for c in cls_nearest_DCD[icls]:
#        f_out.write(' %f' % (c,))
#    f_out.write('\n')
#f_out.write('\n')

f_out.close()

