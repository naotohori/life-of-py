#!/usr/bin/env python

import math

seq = 'UUUUUUUUUU'
f_out = open('out.pdb','w')

n = len(seq)
b = 5.0
l = 5.0

theta = 2 * math.pi / float(n)
r = b / theta

for i in range(n):
    s = seq[i-1]

    x = r * math.cos(theta * i)
    y = 0
    z = r * math.sin(theta * i)
    f_out.write('ATOM  %5i  P   R%s  A %3i    %8.3f%8.3f%8.3f\n' % (3*i+1, s, i+1, x,y,z))

    x = r * math.cos(theta * (i+0.5))
    y = - math.sqrt(l*l - b*b*0.25) 
    z = r * math.sin(theta * (i+0.5))
    f_out.write('ATOM  %5i  S   R%s  A %3i    %8.3f%8.3f%8.3f\n' % (3*i+2, s, i+1, x,y,z))

    x = r * math.cos(theta * (i+0.5))
    y = - math.sqrt(l*l - b*b*0.25) - l
    z = r * math.sin(theta * (i+0.5))
    f_out.write('ATOM  %5i  %sb  R%s  A %3i    %8.3f%8.3f%8.3f\n' % (3*i+3, s, s, i+1, x,y,z))

f_out.close()
