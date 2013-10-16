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
#gnuplotで、全域を描画するために、HをH_pltへ変換する。


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

#dist_bins = [x*5.0 for x in xrange(0,61)]
dist_bins = [x*10.0 for x in xrange(0,31)]

phi = []
theta = [] 
#dist = []
#for i in xrange(11):
#    dist.append([])


### 距離依存の角度に関する統計
## 1: r < 100
## 2: 100 < r < 150
## 3: 150 < r
#phi1=[]
#theta1=[]
#phi2=[]
#theta2=[]
#phi3=[]
#theta3=[]

### 角度依存の距離に関する統計
## N:北半球
## S:南半球
## 0: 極側 (phi依存性なし)
## 1-4: 赤道側は、90度刻みで４つの領域に分ける（番号はGreenwichの裏側から時計回り)
##    1: 135 < phi , phi < -135 (Greenwichの裏側)
##    2: 45 < phi < 135
##    3: -45 < phi < 45 (Greenwichの方角を中心に90度）
##    4: -135 < phi < -45
#dist[1] :N0
#dist[2] :N1
#dist[3] :N2
#dist[4] :N3
#dist[5] :N4
#dist[6] :S1
#dist[7] :S2
#dist[8] :S3
#dist[9] :S4
#dist[10] :S0

############################

for l in file_in:
    if l.find('#') != -1:
        continue
    lsp = l.split()
    r = float(lsp[COL_DIST])
    t = float(lsp[COL_THETA])
    p = float(lsp[COL_PHI])
    #dist[0].append(r)
    theta.append(t)
    phi.append(p)
        
    #if t < theta_bins[1]:
    #    dist[1].append(r)
    #elif t < 90.0:
    #    if p >= 135.0 or p < -135.0:
    #        dist[2].append(r)
    #    elif p > 45.0 and p <= 135.0:
    #        dist[3].append(r)
    #    elif p > -45.0 and p <= 45.0:
    #        dist[4].append(r)
    #    else:
    #        dist[5].append(r)
    #elif t < theta_bins[-2]:
    #    if p >= 135.0 or p < -135.0:
    #        dist[6].append(r)
    #    elif p > 45.0 and p <= 135.0:
    #        dist[7].append(r)
    #    elif p > -45.0 and p <= 45.0:
    #        dist[8].append(r)
    #    else:
    #        dist[9].append(r)
    #else:
    #    dist[10].append(r)
        
################# calc histogram and write to file

############### polar
file_gnu = open("hist_pol.gnudat",'w')

H, theta_edge, phi_edge = histogram2d(theta,phi,bins=[theta_bins,phi_bins])
    #Hd, theta_edge, phi_edge = histogram2d(theta[n],phi[n],bins=[theta_bins,phi_bins],
    #                                       normed=True)
    # normed=Trueにすると変になる。（範囲の広さに応じてnormalizeしてる?)
    
H_plt = convert_array_for3D(H)
    #Hd_plt = convert_array_for3D(Hd)
hsum = H.sum()
    
for i in xrange(len(H_plt)):
    for j in xrange(len(H_plt[i])):
        file_gnu.write('%8.3f %8.3f %15.10e %10.3f\n' % (phi_edge[j],theta_edge[i],H_plt[i,j]/float(hsum),H_plt[i,j]))
    file_gnu.write("\n")
file_gnu.write("\n\n")
        
file_gnu.close()


############### distance
#file_out = open(file_pfx+"_hist_dist.out",'w')
#for d in dist:
#    
#    H, dist_edge = histogram(d,bins=dist_bins)
#    Hd, dist_edge = histogram(d,bins=dist_bins,normed=True)
#    
#    for i,x in enumerate(H):
#        file_out.write('%8.3f %8.6f %10i %8.3f %8.3f\n'
#                        % ((dist_edge[i]+dist_edge[i+1])*0.5, Hd[i], x, dist_edge[i], dist_edge[i+1]))
#    file_out.write("\n\n")
#    
#file_out.close()