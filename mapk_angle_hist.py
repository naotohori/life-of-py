#!/usr/bin/env python
#vim:fileencoding=UTF-8



import sys
import math
from numpy import histogram

if len(sys.argv) != 3:
    print('Usage: SCRIPT [input data] [output prefix]')
    sys.exit(2)
    
file_in = open(sys.argv[1],'r')
file_pfx = sys.argv[2]

COL_THETA = 2 - 1

#theta_bins = [x*10.0 for x in xrange(0,19)] # 10度
theta_bins = [x*5.0 for x in range(0,37)] #  5度
#phi_bins = [x*15.0 for x in xrange(-12,13)] # 15度


theta = [] 
weight = []

for l in file_in:
    if l.find('#') != -1:
        continue
    lsp = l.split()
    t = float(lsp[COL_THETA])
    theta.append(t)
    weight.append(1.0/math.sin(math.radians(t)))
        
        
h, theta_edge = histogram(theta,bins=theta_bins)
hw, theta_edge = histogram(theta,weights=weight, bins=theta_bins)
#hd, theta_edge = histogram(theta,bins=theta_bins,density=True)

hsum = float(sum(h))
hwsum = sum(hw)

file_out = open(file_pfx+"_hist.out",'w')
    
for i,x in enumerate(h):
    ang = (theta_edge[i]+theta_edge[i+1])*0.5
    jcb = math.sin(math.radians(ang))
    file_out.write('%8.3f %10.6f %10.6f %10i\n'
                    % (ang, x/hsum/jcb, hw[i]/hwsum, x))
    
file_out.close()