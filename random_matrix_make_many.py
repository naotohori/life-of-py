#!/usr/bin/env python

import sys
import random
from .matrix_transform_make import matrix_transform_make
from math import pi

if len(sys.argv) != 7:
     print('Usage: % SCRIPT [number of files] [seed] [filename prefix] [x] [y] [z]')
     sys.exit(2)

n = int(sys.argv[1])
random.seed(int(sys.argv[2]))
filename_pre = sys.argv[3]
x = float(sys.argv[4])
y = float(sys.argv[5])
z = float(sys.argv[6])

for i in range(1,n+1) :
     filename = '%s_%04i.mat' % (filename_pre, i)
     a = random.random() * 2.0 * pi
     b = random.random() * 2.0 * pi
     c = random.random() * 2.0 * pi
     matrix_transform_make(a,b,c,x,y,z,filename)
   
