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
#20130702_2.pyよりコピー
#matplotlibは使わない。（補完のしかたが謎）
#gnuplotで、全域を描画するために、HをHpltへ変換する。


import sys
import math
from numpy import histogram, histogram2d
from gnu_data import convert_array_for3D

if len(sys.argv) != 3:
    print 'Usage: SCRIPT [input data] [output prefix]'
    sys.exit(2)
    
file_in = open(sys.argv[1],'r')
file_pfx = sys.argv[2]

COL_DIST = 1 - 1
COL_THETA = 4 - 1
COL_PHI = 5 - 1

phi_bins = [x*10.0 for x in xrange(-18,19)] # 10度
#phi_bins = [x*15.0 for x in xrange(-12,13)] # 15度

DIV_Z = 17  # nの値。-180から180は 2n = 2 x DIV_Z に分割される。
theta_bins = []
theta_bins.append(0.0)          # i=0
for i in xrange(1,DIV_Z):       # i=1,2,3,....,(n-1)
    theta_bins.append(math.degrees(math.acos(1.0-i/float(DIV_Z))))
theta_bins.append(90.0)         # i=n
for i in xrange(1,DIV_Z):       # i=(n+1),(n+2),....,(2n-1)
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

file_gnu = open("hist_pol.gnudat",'w')
##### polar
H, theta_edge, phi_edge = histogram2d(theta,phi,bins=[theta_bins,phi_bins])

Hplt = convert_array_for3D(H)

for i in xrange(len(Hplt)):
    file_gnu.write("\n")
    for j in xrange(len(Hplt[i])):
        file_gnu.write('%f %f %f\n' % (phi_edge[j],theta_edge[i],Hplt[i,j]))

file_gnu.write("\n\n")

##### polar 1: <100
H, theta_edge, phi_edge = histogram2d(theta1,phi1,bins=[theta_bins,phi_bins])

Hplt = convert_array_for3D(H)

file_gnu.write("# hist_pol <100")
for i in xrange(len(Hplt)):
    file_gnu.write("\n")
    for j in xrange(len(Hplt[i])):
        file_gnu.write('%f %f %f\n' % (phi_edge[j],theta_edge[i],Hplt[i,j]))
        
file_gnu.write("\n\n")

##### polar 2: <150
H, theta_edge, phi_edge = histogram2d(theta2,phi2,bins=[theta_bins,phi_bins])

Hplt = convert_array_for3D(H)

file_gnu.write("# hist_pol >100 and <150")
for i in xrange(len(Hplt)):
    file_gnu.write("\n")
    for j in xrange(len(Hplt[i])):
        file_gnu.write('%f %f %f\n' % (phi_edge[j],theta_edge[i],Hplt[i,j]))
        
file_gnu.write("\n\n")
    
##### polar 3: >150
H, theta_edge, phi_edge = histogram2d(theta3,phi3,bins=[theta_bins,phi_bins])

Hplt = convert_array_for3D(H)

file_gnu.write("# hist_pol >150")
for i in xrange(len(Hplt)):
    file_gnu.write("\n")
    for j in xrange(len(Hplt[i])):
        file_gnu.write('%f %f %f\n' % (phi_edge[j],theta_edge[i],Hplt[i,j]))
        
file_gnu.close()

#### distance
dist_bins = [x*5.0 for x in xrange(0,61)]
#H, dist_edge = histogram(dist,bins=dist_bins,normed=True)
H, dist_edge = histogram(dist,bins=dist_bins)
file_out = open(file_pfx+"_dist.out",'w')
for i,x in enumerate(H):
    file_out.write('%8.3f %8.6f %8.3f %8.3f\n'
                    % ((dist_edge[i]+dist_edge[i+1])*0.5, x, dist_edge[i], dist_edge[i+1]))
file_out.close()