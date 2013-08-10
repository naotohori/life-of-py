#!/usr/bin/env python
#vim:fileencoding=UTF-8

#mapk_polar_hist.pyよりコピー

import sys
from numpy import histogram

if len(sys.argv) != 3:
    print 'Usage: SCRIPT [input data] [output prefix]'
    sys.exit(2)
    
file_in = open(sys.argv[1],'r')
file_pfx = sys.argv[2]

COL_DIST = 1 - 1

dist_bins = [x*2.0 for x in xrange(0,81)]

dist = []


for l in file_in:
    if l.find('#') != -1:
        continue
    lsp = l.split()
    r = float(lsp[COL_DIST])
    dist.append(r)
        
############### distance
file_out = open(file_pfx+"_hist_Rg.out",'w')
    
H, dist_edge = histogram(dist,bins=dist_bins)
Hd, dist_edge = histogram(dist,bins=dist_bins,normed=True)
    
for i,x in enumerate(H):
    file_out.write('%8.3f %8.6f %10i %8.3f %8.3f\n'
                    % ((dist_edge[i]+dist_edge[i+1])*0.5, Hd[i], x, dist_edge[i], dist_edge[i+1]))
file_out.write("\n\n")
    
file_out.close()