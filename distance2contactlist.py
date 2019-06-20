#!/usr/bin/env python

import numpy as np
import math
import sys

if len(sys.argv) != 5:
    print('SCRIPT [CG pdb] [pairdist.dat] [Rc] [output prefix]')
    sys.exit(2)

d_cut = float(sys.argv[3])

num_nt = 0
xyz_com = []
for l in open(sys.argv[1],'r'):
    if l[0:4] != 'ATOM':
        continue
    x = float(l[30:38])
    y = float(l[38:46])
    z = float(l[46:54])
    xyz_com.append((x,y,z))

num_nt = len(xyz_com)
print("num_nt: ",num_nt)

dat = np.zeros((num_nt+1, num_nt+1))

ncon = 0
for l in open(sys.argv[2],'r'):
    lsp = l.split()
    if len(lsp) < 3:
        continue
    i = int(lsp[0])
    j = int(lsp[1])
    d = float(lsp[2])
    if abs(i-j) > 2 and d < d_cut:
        dat[int(i), int(j)] = 1
        dat[int(j), int(i)] = 1
        ncon += 1
        
### data
f_out = open(sys.argv[-1] + '.dat','w')
for i in range(1,num_nt+1):
    f_out.write('%i' % (dat[i,1],))
    for j in range(2,num_nt+1):
        f_out.write(' %i' % (dat[i,j],))
    f_out.write('\n')


### For gnuplot
f_out = open(sys.argv[-1] + '.gnudat','w')

f_out.write('%i %i %f\n' % (0,0,dat[1,1]))
for j in range(1,num_nt):
    f_out.write('%i %i %f\n' % (0,j,dat[1,j]))
    f_out.write('%i %i %f\n' % (0,j,dat[1,j+1]))
f_out.write('%i %i %f\n' % (0,num_nt,dat[1,num_nt]))
f_out.write('\n')

for i in range(1,num_nt):
    f_out.write('%i %i %f\n' % (i,0,dat[i,1]))
    for j in range(1,num_nt):
        f_out.write('%i %i %f\n' % (i,j,dat[i,j]))
        f_out.write('%i %i %f\n' % (i,j,dat[i,j+1]))
    f_out.write('%i %i %f\n' % (i,num_nt,dat[i,num_nt]))
    f_out.write('\n')

    f_out.write('%i %i %f\n' % (i,0,dat[i+1,1]))
    for j in range(1,num_nt):
        f_out.write('%i %i %f\n' % (i,j,dat[i+1,j]))
        f_out.write('%i %i %f\n' % (i,j,dat[i+1,j+1]))
    f_out.write('%i %i %f\n' % (i,num_nt,dat[i+1,num_nt]))
    f_out.write('\n')

f_out.write('%i %i %f\n' % (num_nt,0,dat[num_nt,1]))
for j in range(1,num_nt):
    f_out.write('%i %i %f\n' % (num_nt,j,dat[num_nt,j]))
    f_out.write('%i %i %f\n' % (num_nt,j,dat[num_nt,j+1]))
f_out.write('%i %i %f\n' % (num_nt,num_nt,dat[num_nt,num_nt]))
f_out.write('\n')
