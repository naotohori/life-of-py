#!/usr/bin/env python

from numpy import zeros
from math import cos, sin

class mtx_coord_transform():
    def __init__(self):
        self.mtx = zeros((4,4))
        self.mtx[3,3] = 1.0
        
    def translation(self,x,y,z):
        '''並進移動'''
        self.mtx[0,3] = x
        self.mtx[1,3] = y
        self.mtx[2,3] = z
         
    def euler_zxz(self,a,b,c):
        '''Z-X-Z系のオイラー角で回転'''
        self.mtx[0,0] = cos(a)*cos(c) - sin(a)*cos(b)*sin(c)
        self.mtx[1,0] = cos(a)*sin(c) + sin(a)*cos(b)*cos(c)
        self.mtx[2,0] = sin(a)*sin(b)
    
        self.mtx[0,1] = - sin(a)*cos(c) - cos(a)*cos(b)*sin(c)
        self.mtx[1,1] = - sin(a)*sin(c) + cos(a)*cos(b)*cos(c)
        self.mtx[2,1] = cos(a)*sin(b)
    
        self.mtx[0,2] = sin(b)*sin(c)
        self.mtx[1,2] = - sin(b)*cos(c)
        self.mtx[2,2] = cos(b)
    
    def rotate(self, nx, ny, nz, t):
        '''任意の単位ベクトル(nx,ny,nz)を軸として、角度tだけ回転する。'''
        self.mtx[0,0] = nx*nx*(1.0-cos(t)) +    cos(t)
        self.mtx[1,0] = nx*ny*(1.0-cos(t)) + nz*sin(t)
        self.mtx[2,0] = nz*nx*(1.0-cos(t)) - ny*sin(t)
    
        self.mtx[0,1] = nx*ny*(1.0-cos(t)) - nz*sin(t)
        self.mtx[1,1] = ny*ny*(1.0-cos(t)) +    cos(t)
        self.mtx[2,1] = nz*nx*(1.0-cos(t)) + nx*sin(t)
    
        self.mtx[0,2] = nz*nx*(1.0-cos(t)) + ny*sin(t)
        self.mtx[1,2] = ny*nz*(1.0-cos(t)) - nx*sin(t)
        self.mtx[2,2] = nz*nz*(1.0-cos(t)) +    cos(t)
     
         
if __name__ == "__main__" :
    import sys
    
    if not len(sys.argv) in (7,8):
        print ''
        print 'This script makes a homogeneous transformation matrix,'
        print 'angles of which is defined by Z-X-Z Euler angles.'
        print ''
        print 'Usage: % SCRIPT [alpha] [beta] [gamma] [x] [y] [z] [[output]]'
        print ''
        print 'When "output" is specified, the matrix will be written in the file. Otherwise STDOUT is used to display.'
        sys.exit(2)

    a = float(sys.argv[1])
    b = float(sys.argv[2])
    c = float(sys.argv[3])
    x = float(sys.argv[4])
    y = float(sys.argv[5])
    z = float(sys.argv[6])
    
    mtx = mtx_coord_transform()
    mtx.translation(x, y, z)
    mtx.euler_zxz(a, b, c)
    
    if len(sys.argv) == 8:
        file_mat = file(sys.argv[-1],'w')
        file_mat.write('#matrix\n')
        file_mat.write('%15.10f %15.10f %15.10f %15.10f\n' % tuple(mtx.mtx[0]))
        file_mat.write('%15.10f %15.10f %15.10f %15.10f\n' % tuple(mtx.mtx[1]))
        file_mat.write('%15.10f %15.10f %15.10f %15.10f\n' % tuple(mtx.mtx[2]))
        file_mat.write('%15.10f %15.10f %15.10f %15.10f\n' % tuple(mtx.mtx[3]))
        file_mat.write('#a: %f\n#b: %f\n#c: %f\n' % (a,b,c)) 
        file_mat.write('#x: %f\n#y: %f\n#z: %f\n' % (x,y,z)) 
        file_mat.close()
    else:
        sys.stdout.write('#matrix\n')
        sys.stdout.write('%15.10f %15.10f %15.10f %15.10f\n' % tuple(mtx.mtx[0]))
        sys.stdout.write('%15.10f %15.10f %15.10f %15.10f\n' % tuple(mtx.mtx[1]))
        sys.stdout.write('%15.10f %15.10f %15.10f %15.10f\n' % tuple(mtx.mtx[2]))
        sys.stdout.write('%15.10f %15.10f %15.10f %15.10f\n' % tuple(mtx.mtx[3]))
        sys.stdout.write('\n')
        sys.stdout.write('#a: %f\n#b: %f\n#c: %f\n' % (a,b,c)) 
        sys.stdout.write('#x: %f\n#y: %f\n#z: %f\n' % (x,y,z)) 