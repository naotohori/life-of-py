#!/usr/bin/env python

from numpy import zeros
from math import cos, sin
import sys

def matrix_transform_make(a,b,c,x,y,z,filename) :
    transform = zeros((4,4))
    
    transform[0,0] = cos(a)*cos(c) - sin(a)*cos(b)*sin(c)
    transform[1,0] = cos(a)*sin(c) + sin(a)*cos(b)*cos(c)
    transform[2,0] = sin(a)*sin(b)
    transform[3,0] = 0.0
    
    transform[0,1] = - sin(a)*cos(c) - cos(a)*cos(b)*sin(c)
    transform[1,1] = - sin(a)*sin(c) + cos(a)*cos(b)*cos(c)
    transform[2,1] = cos(a)*sin(b)
    transform[3,1] = 0.0
    
    transform[0,2] = sin(b)*sin(c)
    transform[1,2] = - sin(b)*cos(c)
    transform[2,2] = cos(b)
    transform[3,2] = 0.0
     
    transform[0,3] = x
    transform[1,3] = y
    transform[2,3] = z
    transform[3,3] = 1.0
     
    file_mat = file(filename,'w')
    file_mat.write('#matrix\n')
    file_mat.write('%15.10f %15.10f %15.10f %15.10f\n' % tuple(transform[0]))
    file_mat.write('%15.10f %15.10f %15.10f %15.10f\n' % tuple(transform[1]))
    file_mat.write('%15.10f %15.10f %15.10f %15.10f\n' % tuple(transform[2]))
    file_mat.write('%15.10f %15.10f %15.10f %15.10f\n' % tuple(transform[3]))
    file_mat.write('#a: %f\n#b: %f\n#c: %f\n' % (a,b,c)) 
    file_mat.write('#x: %f\n#y: %f\n#z: %f\n' % (x,y,z)) 
    file_mat.close()
         
if __name__ == "__main__" :
    if len(sys.argv) != 8:
        print('')
        print('This script makes a homogeneous transformation matrix,')
        print('angles of which is defined by Z-X-Z Euler angles.')
        print('')
        print('Usage: % SCRIPT [alpha] [beta] [gamma] [x] [y] [z] [output]')
        sys.exit(2)

    a = float(sys.argv[1])
    b = float(sys.argv[2])
    c = float(sys.argv[3])
    x = float(sys.argv[4])
    y = float(sys.argv[5])
    z = float(sys.argv[6])
    filename = sys.argv[-1]
    matrix_transform_make(a,b,c,x,y,z,filename)

