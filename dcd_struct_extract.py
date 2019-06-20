#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@author: Naoto Hori
'''

import sys
import copy
from cafysis.file_io.dcd import DcdFile

if len(sys.argv) < 5:
    print() 
    print('Usage: % SCRIPT [input DCD] [output DCD] [[(mp begin) (mp end)] ....]')
    print() 
    sys.exit(2)
    
# input
f_dcd_in = DcdFile(sys.argv[1])
f_dcd_out = DcdFile(sys.argv[2])

id_pairs = []
n = 1
for arg in sys.argv[3:] :
    if n == 1:
        tp = (int(arg),)
    else :
        id_pairs.append(tp + (int(arg),))
    n *= -1

if n == -1:
    print() 
    print('Usage: % SCRIPT [input DCD] [output DCD] [[(mp begin) (mp end)] ....]')
    sys.exit(2)
    
f_dcd_in.open_to_read()
f_dcd_in.read_header()
header = copy.deepcopy(f_dcd_in.get_header())
nmp_in = header.nmp_real

# Reduce the 'nmp_real' of new header
header.nmp_real = 0
for id_pair in id_pairs:
    header.nmp_real += (id_pair[1] - id_pair[0] + 1)
nmp_out = header.nmp_real
    
f_dcd_out.open_to_write()
f_dcd_out.set_header(header)
f_dcd_out.write_header()
##debug
#print nmp_out

while f_dcd_in.has_more_data() :
    data = f_dcd_in.read_onestep()
    
    data_out = []
    for id_pair in id_pairs :
        if id_pair[1] > nmp_in + 1:
            print('Error: id_pair[1] > nmp_in + 1')
            sys.exit(2)
        data_out.extend(data[id_pair[0]-1 : id_pair[1]])
    ##debug
    #print len(data_out)
    
    f_dcd_out.write_onestep(data_out) 

f_dcd_in.close()
f_dcd_out.close()   
    
