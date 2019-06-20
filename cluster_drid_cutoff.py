#!/usr/bin/env python

import sys
import math
import numpy as np
import scipy.cluster.hierarchy
from cafysis.file_io.drid import DridFile

if len(sys.argv) != 5:
    print('Usage: SCRIPT [DRID file] [prefix] [cutoff] [nskip (to calculate frame id)]')
    sys.exit(2)

drid_filepath = sys.argv[1]
prefix = sys.argv[2]
cutoff_char = sys.argv[3]
cutoff = float(sys.argv[3])
nskip = int(sys.argv[4])

z = []
for l in open(prefix+'.drid.cls.z'):
    lsp = l.split()
    z.append( [int(lsp[1]), int(lsp[2]), float(lsp[3]), int(lsp[4])] )


fcls = scipy.cluster.hierarchy.fcluster(z, cutoff, criterion='distance')   # 300

#print len(fcls)
#print fcls
f_out = open(prefix+'.drid.cls_%s' % cutoff_char, 'w')
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

#fs = []
#for i in range(ncls):
#    fs.append(open(prefix+".drid.cls.%i" % (i+1,),'w'))
#
#f_traj = open("%s.e_rmsd_ddrid_pca_hbnn" % (prefix,), 'r')
#
#for i in range(len(fcls)):
#    l = f_traj.readline()
#    icls = fcls[i]
#    fs[icls-1].write(l)
#    for iskip in range(nskip-1):
#        f_traj.readline()
#
#f_traj.close()
#for f in fs:
#    f.close()

f_out = open(prefix+".drid.cls_%s.summary" % cutoff_char, 'w')
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

drid = DridFile(drid_filepath)
drid.open_to_read()
drid.read_header()

cls_centroids = []
nc = 3*drid.get_header().n_centroids()
for icls in range(ncls):
    cls_centroids.append( np.zeros((nc,)) )

k = 0
while drid.has_more_data():
    cls_centroids[fcls[k]-1] += drid.read_onestep()
    k+=1

for icls in range(ncls):
    #cls_centroids[icls] /= float(k)
    cls_centroids[icls] /= float( fcls.tolist().count(icls+1) )

cls_nearest_DRID = []
cls_nearest_node = []
cls_nearest_dDRID = []
cls_average_dDRID = []
cls_num_node = []
for icls in range(ncls):
    cls_nearest_DRID.append( np.zeros((nc,)) )
    cls_nearest_node.append(-1)
    cls_nearest_dDRID.append(9999999.)
    cls_average_dDRID.append(0.0)
    cls_num_node.append(0)

drid.rewind()
k = 0
while drid.has_more_data():
    data = drid.read_onestep()
    icls = fcls[k] - 1
    d = math.sqrt(np.square(cls_centroids[icls][:] - data[:]).sum())
    if d < cls_nearest_dDRID[icls]:
        cls_nearest_dDRID[icls] = d
        cls_nearest_DRID[icls][:] = data[:]
        #cls_nearest_node[icls] = k + 1
        cls_nearest_node[icls] = k
    cls_average_dDRID[icls] += d
    cls_num_node[icls] += 1
    k += 1

# Check this code correctly working
for icls in range(ncls):
    if cls_num_node[icls] != fcls.tolist().count(icls+1):
        print('Error: cls_num_node[icls] != fcls.count(icls+1)')
        sys.exit(2)

for icls in range(ncls):
    cls_average_dDRID[icls] /= float(cls_num_node[icls])

f_out = open('%s.drid.cls_%s.centroids' % (prefix, cutoff_char),'w')
f_out.write('#clsID node dDRID dDRID/sqrt(nc) <dDRID> <dDRID>/sqrt(nc)\n')
for icls in range(ncls):
    f_out.write('%i %i %f %f %f %f\n' % (icls+1, cls_nearest_node[icls], 
            cls_nearest_dDRID[icls], cls_nearest_dDRID[icls]/math.sqrt(nc),
            cls_average_dDRID[icls], cls_average_dDRID[icls]/math.sqrt(nc)) )
f_out.write('\n')

f_out.write('#centroid DRID\n')
for icls in range(ncls):
    f_out.write('%i' % (icls+1,))
    for c in cls_centroids[icls]:
        f_out.write(' %f' % (c,))
    f_out.write('\n')
f_out.write('\n')

f_out.write('#nearest DRID\n')
for icls in range(ncls):
    f_out.write('%i' % (icls+1,))
    for c in cls_nearest_DRID[icls]:
        f_out.write(' %f' % (c,))
    f_out.write('\n')
f_out.write('\n')

f_out.close()

