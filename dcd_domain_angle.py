#!/usr/bin/env python

'''
See pdb_domain_angle.py for detail and test.

To calculate domain angle (gamma) in Supplementary Figure 7
in Denesyuk & Thirumalai, Nat Chem (2015)
'''

'''
Domain A P5-P4-P6a:   52 - 124
Domain B P7-P3-P8:    42 -  47
                     174 - 178
                     130 - 166
'''

import sys
import numpy as np
import math
from lop.file_io.dcd import DcdFile

if len(sys.argv) != 3:
    print ('Usage: SCRIPT [DCD file] [angle file]')
    sys.exit(2)

f_out = open(sys.argv[-1],'w')

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


dcd = DcdFile(sys.argv[1])
dcd.open_to_read()
dcd.read_header()

while dcd.has_more_data():
    data = dcd.read_onestep()

    ''' Calculate principal axis of domain A '''
    xyzs = []
    for imp in mp_domA:
        xyzs.append( data[imp-1] )

    xyzs = np.array(xyzs)
    com = np.sum(xyzs, axis=0) / float(xyzs.shape[0])
    xyzs = xyzs - com

    cov = np.cov(xyzs, rowvar=False)
    _, _, vA =  np.linalg.svd(cov)

    r0 = np.dot(data[ mp_domA_0-1 ], vA[0])
    r1 = np.dot(data[ mp_domA_1-1 ], vA[0])
    if r1 < r0:
        vA[0] *= -1.0


    ''' Calculate principal axis of domain B '''
    xyzs = []
    for imp in mp_domB:
        xyzs.append( data[imp-1] )

    xyzs = np.array(xyzs)
    com = np.sum(xyzs, axis=0) / float(xyzs.shape[0])
    xyzs = xyzs - com

    cov = np.cov(xyzs, rowvar=False)
    _, _, vB =  np.linalg.svd(cov)

    r0 = np.dot(data[ mp_domB_0-1 ], vB[0])
    r1 = np.dot(data[ mp_domB_1-1 ], vB[0])
    if r1 < r0:
        vB[0] *= -1.0


    ''' Angle between domain A and B '''
    #cos_theta = np.dot(vA[0],vB[0]) / sqrt(np.dot(vA[0],vA[0]) * np.dot(vB[0],vB[0]))
    theta = math.acos( np.dot(vA[0],vB[0]) / math.sqrt(np.dot(vA[0],vA[0]) * np.dot(vB[0],vB[0])) )
    theta =  180.0 * theta / math.pi

    f_out.write('%6.2f\n' %  theta)
dcd.close()
