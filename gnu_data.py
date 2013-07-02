#!/usr/bin/env python
#vim:fileencoding=UTF-8

def convert_array_for3D(H):
    '''gnuplotで、全域を描画するために、HをHpltへ変換する。'''
    
    from numpy import zeros
    
    n = len(H)
    m = len(H[0])
    Hplt = zeros((n+1,m+1))
    
    Hplt[0,0] = H[0,0]
    for i in xrange(1, n):
        Hplt[i,0] = (H[i-1,0]+H[i,0]) * 0.5
    Hplt[n,0] = H[n-1,0]
    
    for j in xrange(1,m):
        Hplt[0,j] = (H[0,j-1]+H[0,j]) * 0.5
        for i in xrange(1,n):
            Hplt[i,j] = (H[i-1,j-1]+H[i-1,j]+H[i,j-1]+H[i,j]) * 0.25
        Hplt[n,j] = (H[n-1,j-1]+H[n-1,j]) * 0.5
    
    Hplt[0,m] = H[0,m-1]
    for i in xrange(1,n):
        Hplt[i,m] = (H[i-1,m-1]+H[i,m-1]) * 0.5
    Hplt[n,m] = H[n-1,m-1]
                
    return Hplt