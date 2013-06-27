#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2013/03/12
@author: Naoto Hori
'''

import sys
#import math
from cafysis.file_io.dcd import DcdFile

if not len(sys.argv) in (5,):
    print 'Currently this script is available for a system which contains only one RNA molecule.'
    print 'Usage: % SCRIPT [input DCD] [first particle ID for calc] [gap] [output]'
    sys.exit(2)
    
offset = int(sys.argv[2])
gap = int(sys.argv[3])
filepath_out = sys.argv[4]

dcd = DcdFile(sys.argv[1])
dcd.open_to_read()
dcd.read_header()

nmp = dcd.get_header().nmp_real

unit_len_sq = 0.0
sum_cos_theta = [0.0] * nmp
num_n = [0] * nmp

first = True
while dcd.has_more_data() :
    data = dcd.read_onestep()
    if first:
        first = False
        unit_len_sq = ( (data[offset-1][0] - data[offset-1+gap][0]) ** 2
                       +(data[offset-1][1] - data[offset-1+gap][1]) ** 2
                       +(data[offset-1][2] - data[offset-1+gap][2]) ** 2 )
        
    for i in xrange(offset, nmp-2, gap) :
        for j in xrange(i+gap, nmp-1, gap) :
            # i,j
            ij = (  data[i][0] * data[j][0]
                  + data[i][1] * data[j][1]
                  + data[i][2] * data[j][2] ) 
            # i-1,j
            i_1_j = (  data[i-1][0] * data[j][0]
                     + data[i-1][1] * data[j][1]
                     + data[i-1][2] * data[j][2] ) 
            # i,j-1
            i_j_1 = (  data[i][0] * data[j-1][0]
                     + data[i][1] * data[j-1][1]
                     + data[i][2] * data[j-1][2] ) 
            # i-1,j-1
            i_1_j_1 = (  data[i-1][0] * data[j-1][0]
                       + data[i-1][1] * data[j-1][1]
                       + data[i-1][2] * data[j-1][2] ) 
            cos_theta = ( ij - i_1_j - i_j_1 + i_1_j_1 ) / unit_len_sq
            n = (j-i) / gap 
            sum_cos_theta[n] += cos_theta
            num_n[n] += 1
            #print n, cos_theta
            
    for i in xrange(nmp/gap-2):
    #for i in xrange(10):
        if num_n[i] == 0:
            print i, 1.0
        else:
            print i, sum_cos_theta[i]/float(num_n[i])
            #print i, (- i * math.sqrt(unit_len_sq) / math.log(sum_cos_theta[i]/float(num_n[i])))
    print ''
    print ''
            