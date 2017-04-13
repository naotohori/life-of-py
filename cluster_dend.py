#!/usr/bin/env python

import sys
import scipy.cluster.hierarchy

if len(sys.argv) not in (3,4):
    print 'Usage: SCRIPT [prefix] [p (number of leaves)]'
    print 'Usage: SCRIPT [prefix] [p (number of leaves)] [["NOLABEL"]]'
    sys.exit(2)

if len(sys.argv) == 3:
    flg_no_label = False
else:
    if sys.argv[3] == "NOLABEL":
        flg_no_label = True
    else:
        print 'Usage: SCRIPT [prefix] [p (number of leaves)]'
        print 'Usage: SCRIPT [prefix] [p (number of leaves)] [["NOLABEL"]]'
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
                                no_labels=flg_no_label,
                                show_leaf_counts=True)
plt.savefig(prefix+'.drid.cls.dend_%i.ps' % p)
plt.figure()
plt.clf()


