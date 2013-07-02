#!/usr/bin/env python
#vim:fileencoding=UTF-8

import sys
from numpy import histogram, histogram2d, zeros
import math
import matplotlib.pyplot as plt
from scipy.constants.constants import pi

if len(sys.argv) != 4:
    print 'Usage: SCRIPT [input data] [output prefix] [output file]'
    sys.exit(2)
    
file_in = open(sys.argv[1],'r')
file_pfx = sys.argv[2]
file_out = open(sys.argv[3],'w')

COL_DIST = 1 - 1
COL_THETA = 4 - 1
COL_PHI = 5 - 1

dist = []
phi = []
theta = [] 

phi1=[]
theta1=[]
phi2=[]
theta2=[]
phi3=[]
theta3=[]

for l in file_in:
    if l.find('#') != -1:
        continue
    lsp = l.split()
    r = float(lsp[COL_DIST])
    t = float(lsp[COL_THETA])
    p = float(lsp[COL_PHI])
    dist.append(r)
    theta.append(t)
    phi.append(p)
    if r < 100.0:
        theta1.append(t)
        phi1.append(p)
    elif r < 150.0:
        theta2.append(t)
        phi2.append(p)
    else:
        theta3.append(t)
        phi3.append(p)
    
theta_bins = [x*10.0 for x in xrange(0,19)]
phi_bins = [x*10.0 for x in xrange(-18,19)]
#theta_bins = [x*15.0 for x in xrange(0,13)]
#phi_bins = [x*15.0 for x in xrange(-12,13)]

weight = []
EPSI = 0.001
for t in theta:
    radt = math.radians(t)
#    weight.append(1.0/math.sin(radt))
    if radt<EPSI or (pi-radt)<EPSI:
        weight.append(1.0/math.sin(EPSI))
    else:
        weight.append(1.0/math.sin(radt))

H, theta_edge, phi_edge = histogram2d(theta,phi,bins=[theta_bins,phi_bins],
                                      normed=False)

Hn, theta_edge, phi_edge = histogram2d(theta,phi,bins=[theta_bins,phi_bins],
                                      normed=True)

Hw, theta_edge, phi_edge = histogram2d(theta,phi,bins=[theta_bins,phi_bins],
                                      normed=False, weights=weight)

Hnw, theta_edge, phi_edge = histogram2d(theta,phi,bins=[theta_bins,phi_bins],
                                      normed=True, weights=weight)

Hj = zeros((18,36))
Hnj = zeros((18,36))
for t in xrange(18):
    t_deg = t*10.0 + 5.0
    for p in xrange(36):
        Hj[t,p] = H[t,p] / math.sin(math.radians(t_deg))
        Hnj[t,p] = Hn[t,p] / math.sin(math.radians(t_deg))
#Hj = zeros((12,24))
#Hnj = zeros((12,24))
#for t in xrange(12):
#    t_deg = t*15.0 + 7.5
#    for p in xrange(24):
#        Hj[t,p] = H[t,p] / math.sin(math.radians(t_deg))
#        Hnj[t,p] = Hn[t,p] / math.sin(math.radians(t_deg))
        


#H1, theta_edge, phi_edge = histogram2d(theta1,phi1,bins=[theta_bins,phi_bins],
#                                      normed=True)
#H2, theta_edge, phi_edge = histogram2d(theta2,phi2,bins=[theta_bins,phi_bins],
#                                      normed=True)
#H3, theta_edge, phi_edge = histogram2d(theta3,phi3,bins=[theta_bins,phi_bins],
#                                      normed=True)

extent=[phi_edge[-1],phi_edge[0],theta_edge[-1],theta_edge[0]]

plt.imshow(H,extent=extent, interpolation='nearest')
plt.colorbar()
plt.savefig(file_pfx+'H.png')

plt.clf()
plt.imshow(Hn,extent=extent, interpolation='nearest')
plt.colorbar()
plt.savefig(file_pfx+'Hn.png')

plt.clf()
#plt.imshow(Hw,extent=extent, interpolation='nearest')
plt.imshow(Hw,extent=extent, interpolation='nearest',vmax=300)
plt.colorbar()
plt.savefig(file_pfx+'Hw.png')

plt.clf()
plt.imshow(Hnw,extent=extent, interpolation='nearest')
plt.colorbar()
plt.savefig(file_pfx+'Hnw.png')

plt.clf()
plt.imshow(Hj,extent=extent, interpolation='nearest')
plt.colorbar()
plt.savefig(file_pfx+'Hj.png')

plt.clf()
plt.imshow(Hnj,extent=extent, interpolation='nearest')
plt.colorbar()
plt.savefig(file_pfx+'Hnj.png')

#plt.imshow(H,extent=extent, interpolation='nearest')
#plt.colorbar()
#plt.savefig(file_pfx+'hist.png')
#
#plt.imshow(H1,extent=extent, interpolation='nearest')
#plt.savefig(file_pfx+'hist1.png')
#
#plt.imshow(H2,extent=extent, interpolation='nearest')
#plt.savefig(file_pfx+'hist2.png')
#
#plt.imshow(H3,extent=extent, interpolation='nearest')
#plt.savefig(file_pfx+'hist3.png')


#dist_bins = [x*10.0 for x in xrange(0,31)]
#H, dist_edge = histogram(dist,bins=dist_bins,normed=True)
#for i,x in enumerate(H):
#    file_out.write('%8.3f %8.6f %8.3f %8.3f\n'
#                    % ((dist_edge[i]+dist_edge[i+1])*0.5, x, dist_edge[i], dist_edge[i+1]))
#
file_out.close()