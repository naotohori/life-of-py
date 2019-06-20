#!/usr/bin/env python
#vim:fileencoding=UTF-8

#20130630_1.pyよりコピー
#Z軸方向を等分割(水平方向の輪切り）することで、ヤコビアンを不要にする。
#
# theta（天頂角）= -180 から 180 度の範囲を、2nに分割するとする。
# 
# i = 1, n （北半球）の範囲では、
#      1 - (i/n) = cos(theta_i)
# ===>  theta_i = arccos(1 - i/n)
# 
# i=1の天頂角をなるべく小さくとりたいと思うところだが、
# するとnが大きくなり，赤道付近での角度が非常に小さくなってしまうので、
# 兼ね合いで調節する必要がある。
#
# theta = 10 deg のとき、nは約65.8
# theta = 15 deg のとき、nは約29.3
# theta = 20 deg のとき、nは約16.6
#
# n=30にすると、theta_1 = 14.8 deg、 theta_29 = 88.1 deg
# n=17にすると、theta_1 = 19.7 deg、 theta_16 = 86.6 deg
#


import sys
from numpy import histogram, histogram2d, zeros
import math
import matplotlib.pyplot as plt
from scipy.constants.constants import pi

if len(sys.argv) != 4:
    print('Usage: SCRIPT [input data] [output prefix] [output file]')
    sys.exit(2)
    
file_in = open(sys.argv[1],'r')
file_pfx = sys.argv[2]
file_out = open(sys.argv[3],'w')

COL_DIST = 1 - 1
COL_THETA = 4 - 1
COL_PHI = 5 - 1

phi_bins = [x*10.0 for x in range(-18,19)]
#phi_bins = [x*15.0 for x in xrange(-12,13)]

DIV_Z = 17  # nの値。-180から180は 2n = 2 x DIV_Z に分割される。
theta_bins = []
theta_bins.append(0.0)          # i=0
for i in range(1,DIV_Z):       # i=1,2,3,....,(n-1)
    theta_bins.append(math.degrees(math.acos(1.0-i/float(DIV_Z))))
theta_bins.append(90.0)         # i=n
for i in range(1,DIV_Z):       # i=(n+1),(n+2),....,(2n-1)
    j = DIV_Z - i
    theta_bins.append(180.0-theta_bins[j])
theta_bins.append(180.0)        # i=2n

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

H, theta_edge, phi_edge = histogram2d(theta,phi,bins=[theta_bins,phi_bins])

#for i in len(H):
#    file_out.write('\n')
#    file_out.write('%f %f %f\n' % (phi_edge[0],theta_edge[i],H[i,0]))
#    for j in len(H[0]):
#        file_out.write('%f %f %i\n'
#                       % ( (phi_edge[j]+phi_edge[j+1])*0.5,
#                           (theta_edge[i]+theta_edge[i+1])*0.5,
#                           H[i,j],))
        
n = len(H)
m = len(H[0])
Hplt = zeros((n+1,m+1))

Hplt[0,0] = H[0,0]
for i in range(1, n):
    Hplt[i,0] = (H[i-1,0]+H[i,0]) * 0.5
Hplt[n,0] = H[n-1,0]

for j in range(1,m):
    Hplt[0,j] = (H[0,j-1]+H[0,j]) * 0.5
    for i in range(1,n):
        Hplt[i,j] = (H[i-1,j-1]+H[i-1,j]+H[i,j-1]+H[i,j]) * 0.25
    Hplt[n,j] = (H[n-1,j-1]+H[n-1,j]) * 0.5

Hplt[0,m] = H[0,m-1]
for i in range(1,n):
    Hplt[i,m] = (H[i-1,m-1]+H[i,m-1]) * 0.5
Hplt[n,m] = H[n-1,m-1]
    
for i in range(len(Hplt)):
    file_out.write("\n")
    for j in range(len(Hplt[i])):
        file_out.write('%f %f %f\n' % (phi_edge[j],theta_edge[i],Hplt[i,j]))
    
        
#file_out.write('#\n')
#file_out.write('#\n')
#file_out.write('#\n')
         

plt.hist2d(phi, theta, bins=[phi_bins,theta_bins])
plt.axis()
plt.colorbar()
plt.savefig(file_pfx+'H.svg')
plt.savefig(file_pfx+'H.png')


#dist_bins = [x*5.0 for x in xrange(0,61)]
#H, dist_edge = histogram(dist,bins=dist_bins,normed=True)
#for i,x in enumerate(H):
#    file_out.write('%8.3f %8.6f %8.3f %8.3f\n'
#                    % ((dist_edge[i]+dist_edge[i+1])*0.5, x, dist_edge[i], dist_edge[i+1]))

file_out.close()