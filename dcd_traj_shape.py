#!/usr/bin/env python
#vim:fileencoding=UTF-8
'''
Created on 2016/11/08
@author: Naoto Hori

Calculate shape parameters, D and S, from coordinates in PDB

T: Inertia tensor
   From T, Rg is also calculated as Rg = sqrt(tr(T))

D: Sphericity  (0 <= D <= 1)
   D = 0 --> perfect sphere
   D > 0 --> anisotropic
S: Spheroidal shape  (-1/4 <= S <= 2)
   S < 0 --> oblate
   S > 0 --> prolate

For perfect spheres, D = S = 0.

Ref:
Dima and Thirumalai, J. Phys. Chem. B 2004, 108: 6564-6570
'''

import sys
from cafysis.file_io.dcd import DcdFile
import numpy as np
import math

if len(sys.argv) != 2:
    print('Usage: SCRIPT [dcd]')
    print('')
    print('Returns: %5i %6.3f %6.3f %6.3f % (N, Rg, D, S)')
    print('    N  = number of particles (atoms)')
    print('    Rg = radius of gyration')
    print('    D  = Sphericity (0 <= D <= 1)')
    print('          D = 0 --> perfect sphere')
    print('          D > 0 --> anisotropic')
    print('    S  = Spheroidal shape (-1/4 <= S <= 2)')
    print('          S < 0 --> oblate')
    print('          S > 0 --> prolate')
    sys.exit(2)



dcd = DcdFile(sys.argv[1])
dcd.open_to_read()
dcd.read_header()


'''
xyzs for N-atom PDB
 [ [ 1x, 1y, 1z],
   [ 2x, 2y, 2z],
   [ 3x, 3y, 3z],
   [ 4x,  .
          .
          .
   [ Nx, Ny, Nz] ]
'''

N = dcd._header.nmp_real
fN = float(N)

def calc_shape(xyzs):
    # Inertia tensor
    T = np.zeros((3,3)) 
    
    for i in range(N):
        for j in range(N):
            for a in range(3):
                for b in range(3):
                    T[a,b] += (xyzs[i][a] - xyzs[j][a]) * (xyzs[i][b] - xyzs[j][b])
    
    T = T / (2.0 * fN**2)
    
    w, v = np.linalg.eig( T )
    
    trT = sum(w)
    #print 'w', w
    #print 'v',v
    #print 'trT = ', trT
    Rg = math.sqrt(trT)
    
    w_avg = trT / 3.0
    
    D = 1.5 * ( (w[0]-w_avg)**2 + (w[1]-w_avg)**2 + (w[2]-w_avg)**2 ) / (trT**2)
    
    S = 27.0 * ( (w[0]-w_avg) * (w[1]-w_avg) * (w[2]-w_avg) ) / (trT**3)
    
    return Rg, D, S

print("#N: ", N)

while dcd.has_more_data():
    data = dcd.read_onestep()
    Rg, D, S = calc_shape(data)

    print(('%6.3f %6.3f %6.3f' % (Rg, D, S)))
