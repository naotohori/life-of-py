#!/usr/bin/env python
# 2011/03/17 coded by Naoto HORI

import copy
import sys
from cafysis.file_io.pdb import PdbFile
from scipy import linalg
from numpy import zeros, float32
import glob

file_out = open('pca_pdb.out', 'w')

data = []
nmp = 0

for f in glob.glob('./*pca.pdb'):

    xyzs = []
    chains = PdbFile(f,'r').read_all_and_close()
    for c in chains:
        for r in c.residues:
            xyzs.append( r.atoms[0].xyz.get_as_list() )

    data.append(xyzs)

    if nmp == 0:
        nmp = len(c.residues)
    else:
        if len(c.residues) != nmp:
            print('len(c.residues) != nmp')
            sys.exit(2)

num_model = len(data)

print('nmp: {:d}'.format(nmp))
print('num_model: {:d}'.format(num_model))

mean = []
mean = copy.deepcopy(data[0])
for d in data[1:]:
    for i in range(nmp):
        for j in range(3):
            mean[i][j] += d[i][j]

file_out.write('#mean\n')
for i in range(nmp):
    for j in range(3):
        mean[i][j] /= float(num_model)
        file_out.write('%20e\n' % mean[i][j]) 


for d in data:
    for i in range(nmp):
        for j in range(3):
            d[i][j] -= mean[i][j]

n = nmp*3
covariance = zeros((n, n), dtype=float32)
variance = zeros((n,), dtype=float32)

for d in data:
    
#    for i_mp in range(nmp) :
#        for i_xyz in range(3) :
#            i = 3 * i_mp + i_xyz
#            for j_mp in range(i_mp+1) :
#                for j_xyz in range(3):
#                    j = 3 * j_mp + j_xyz
#                    covariance[i,j] += d[i_mp][i_xyz] * d[j_mp][j_xyz]
#            variance[i] += d[i_mp][i_xyz]
    for i in range(n):
        i_mp = i // 3
        i_xyz = i % 3
        
        variance[i] += d[i_mp][i_xyz]

        for j in range(i+1):
            j_mp = i // 3
            j_xyz = i % 3

            covariance[i,j] += d[i_mp][i_xyz] * d[j_mp][j_xyz]
            
a = zeros((n,n), dtype=float32)
fn = float(num_model - 1)
for i in range(n):
    for j in range(i+1) :
        a[i,j] = covariance[i,j] / fn - variance[i] * variance[j] / (fn **2)
    for j in range(i+1, n) :
        a[i,j] = covariance[j,i] / fn - variance[i] * variance[j] / (fn **2)

for i in range(nmp*3):
    for j in range(i):
        if a[i,j] != a[j,i]:
            print(i,j)

(w, v) = linalg.eigh(a)
#(w, v) = linalg.eig(a)

#print(w) 
#print(v)

idx = w.argsort()[::-1]   
w = w[idx]
v = v[:,idx]

file_out.write('#value\n')
for i in range(nmp*3) :
    file_out.write('%20e\n' % w[i])

for i in range(nmp*3) :
    file_out.write('#vector%10i\n' % i)
    for j in range(nmp*3) :
        file_out.write('%20e\n' % v[j,i])

file_out.close()
