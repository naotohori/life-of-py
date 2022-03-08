#!/usr/bin/env python

'''
To calculate domain angle (gamma) in Supplementary Figure 7
in Denesyuk & Thirumalai, Nat Chem (2015)

To see the angles,
$ vmd -e azo_domain_angle.vmd
[Tk console] $ source azo_domain_angle.tcl
'''

'''
Domain A P5-P4-P6a:   52 - 124
Domain B P7-P3-P8:    42 -  47
                     174 - 178
                     130 - 166
'''

import sys
from lop.file_io.pdb import PdbFile
import numpy as np
import math

#if len(sys.argv) != 3:
#    print ('Usage: SCRIPT [pdb file] [angle file]')
#    sys.exit(2)

pdbfilepath = './azo.rna.pdb'

mp_domA = []  # mp starts with 1
for nt12 in range(52, 124+1):
    nt = nt12 - 11
    mp_domA.append( 3 * (nt-1) )  # Use phosphate 

mp_domA_0 = 3 * (110 - 11)
mp_domA_1 = 3 * ( 72 - 11)

mp_domB = []
for nt12 in list(range(42, 47+1))+list(range(174,178+1))+list(range(130,166+1)):
    nt = nt12 - 11
    mp_domB.append( 3 * (nt-1) )  # Use phosphate

mp_domB_0 = 3 * (155 - 11)
mp_domB_1 = 3 * (130 - 11)

pdb = PdbFile(pdbfilepath)
pdb.open_to_read()
chains = pdb.read_all()
c = chains[0]
pdb.close()

'''
Domain A
'''
xyzs = []
for imp in mp_domA:
    xyzs.append( c.get_atom( imp-1 ).xyz.get_as_list() )

xyzs = np.array(xyzs)
com = np.sum(xyzs, axis=0) / float(xyzs.shape[0])
#print com
xyzs = xyzs - com

cov = np.cov(xyzs, rowvar=False)
#cov = np.cov(xyzs.T, rowvar=True)  # This gives same result
#print cov

u,s,vA =  np.linalg.svd(cov)

#print u
#print s
#print v

r0 = np.dot(c.get_atom( mp_domA_0-1 ).xyz.get_as_ndarray(), vA[0])
r1 = np.dot(c.get_atom( mp_domA_1-1 ).xyz.get_as_ndarray(), vA[0])

if r1 < r0:
    vA[0] *= -1.0

print('Principal vector of domain A', vA[0])


'''
Domain B
'''
xyzs = []
for imp in mp_domB:
    xyzs.append( c.get_atom( imp-1 ).xyz.get_as_list() )

xyzs = np.array(xyzs)
com = np.sum(xyzs, axis=0) / float(xyzs.shape[0])
#print com
xyzs = xyzs - com

cov = np.cov(xyzs, rowvar=False)
#cov = np.cov(xyzs.T, rowvar=True)  # This gives same result
#print cov

u,s,vB =  np.linalg.svd(cov)

#print u
#print s
#print v
#

r0 = np.dot(c.get_atom( mp_domB_0-1 ).xyz.get_as_ndarray(), vB[0])
r1 = np.dot(c.get_atom( mp_domB_1-1 ).xyz.get_as_ndarray(), vB[0])

if r1 < r0:
    vB[0] *= -1.0

print('Principal vector of domain B', vB[0])


''' Angle between domain A and B '''
#cos_theta = np.dot(vA[0],vB[0]) / math.sqrt(np.dot(vA[0],vA[0]) * np.dot(vB[0],vB[0]))
theta = math.acos( np.dot(vA[0],vB[0]) / math.sqrt(np.dot(vA[0],vA[0]) * np.dot(vB[0],vB[0])) )
print(180.0 * theta / math.pi)
