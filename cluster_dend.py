#!/usr/bin/env python

import sys
import scipy.cluster.hierarchy

if len(sys.argv) != 3:
    print 'Usage: SCRIPT [prefix] [p (number of leaves)]'
    sys.exit(2)

prefix = sys.argv[1]
p = int(sys.argv[2])

z = []
for l in open(prefix+'.drid.cls.z'):
    lsp = l.split()
    z.append( [int(lsp[1]), int(lsp[2]), float(lsp[3]), int(lsp[4])] )


import matplotlib.pyplot as plt
plt.figure()
plt.clf()
scipy.cluster.hierarchy.dendrogram(z, p=p, truncate_mode='lastp',
                                orientation='left', color_threshold=0,
                                show_leaf_counts=True)
plt.savefig(prefix+'.drid.cls.dend_%i.ps' % p)
plt.figure()
plt.clf()


