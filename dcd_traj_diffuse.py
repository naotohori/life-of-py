#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2014/11/03
@author: Naoto Hori
'''
import sys
import math
from cafysis.file_io.dcd import DcdFile

if len(sys.argv) != 4:
    print ('\n Usage: SCRIPT [input DCD file] [interval] [output diffuse file]\n')
    sys.exit(2)
    

dcd = DcdFile(sys.argv[1])

interval = int(sys.argv[2])

nmp = dcd.get_header().nmp_real

f_out = open(sys.argv[-1],'w')

# count the number of steps
dcd.open_to_read()
dcd.read_header()
nstep = 0
com = []
while dcd.has_more_data():
    data = dcd.read_onestep()
    nstep += 1
    s = [0.0, 0.0, 0.0]
    for xyz in data:
        s[0] += xyz[0]
        s[1] += xyz[1]
        s[2] += xyz[2]
    s[0] /= len(data)
    s[1] /= len(data)
    s[2] /= len(data)
    com.append(s)
dcd.close()

if nstep != len(com):
    print 'Error: nstep != len(com)'
    sys.exit(2)
Nmax = int(nstep / 2)
print 'nstep: ',nstep
print 'Nmax: ', Nmax


msd = []
n_add = []  ## just for check (delete later)
msd.append(99999999.9) # index 0 is dummy
n_add.append(-1)       # index 0 is dummy
for dt in range(1, Nmax+1):
    msd.append(0.0)
    n_add.append(0)


for i_orig in range(0, Nmax+1):

    xyz_0 = com[i_orig]
    
    #print 'i_orig: ',i_orig

    for dt in range(1, Nmax+1):
        
        xyz_t = com[i_orig + dt]

        msd[dt] += ((xyz_0[0] - xyz_t[0]) ** 2
                   +(xyz_0[1] - xyz_t[1]) ** 2
                   +(xyz_0[2] - xyz_t[2]) ** 2)
        n_add[dt] += 1

for dt in range(1, Nmax+1):
    f_out.write('%i %f\n' % (dt, msd[dt]/float(n_add[dt])) )
print n_add
