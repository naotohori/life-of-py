#!/usr/bin/env python
'''
@author: Naoto Hori
'''

import sys
from lop.file_io.ninfo import NinfoFile
from lop.elements.ninfo import NinfoSet

if len(sys.argv) != 6 :
    print ('Usage: SCRIPT [ninfo file] [edc file] [dist(aln) file] [factor] [output ninfo file]')
    sys.exit(2)
    
f_ninfo = NinfoFile(sys.argv[1])
f_edc = open(sys.argv[2], 'r')
f_dist = open(sys.argv[3],'r')
factor = float(sys.argv[4])
f_out_ninfo = NinfoFile(sys.argv[-1])

f_ninfo.open_to_read()
ninfo = NinfoSet()
f_ninfo.read_all(ninfo)
f_ninfo.close()

f_out_ninfo.open_to_write()
f_out_ninfo._file.write('#0: '+sys.argv[0]+'\n')
f_out_ninfo._file.write('#1: '+sys.argv[1]+'\n')
f_out_ninfo._file.write('#2: '+sys.argv[2]+'\n')
f_out_ninfo._file.write('#3: '+sys.argv[3]+'\n')
f_out_ninfo._file.write('#4: '+sys.argv[4]+'\n')
f_out_ninfo._file.write('#5: '+sys.argv[5]+'\n')

aln = {}
# Read from dist file
'''
# Jenner        Dunkle
# res    mp     res   mp    dist
     1      1     4     11   1.55
     1      2     4     12   8.02
     2      3     5     13   0.93
'''
for line in f_dist :
    if line.find('#') != -1 :
        continue 
    if line == '\n' :
        continue
    linesp = line.split()
    #     Jenner             Dunkle
    imp_Jenner = int(linesp[1])
    imp_Dunkle = int(linesp[3])
    if imp_Jenner in aln :
        aln[imp_Jenner] = aln[imp_Jenner] + (imp_Dunkle,)
    else :
        aln[imp_Jenner] = (imp_Dunkle,)

energys = {}
for line in f_edc :
    if line.find('#') != -1:
        continue
    if line == '\n' :
        continue
#                     0    1   2   3  4   5   6  7  8   9   10  11 12  13 14 15  
#        f_out.write('%3i %6i %4s %6i %3s | %3i %6i %4s %6i %3s | %6i %6i | %f\n' % info)
    linesp = line.split()
    imp1_orig = int(linesp[3])
    imp2_orig = int(linesp[9])
    if imp1_orig in aln:
        imp1_tp = aln[imp1_orig]
    else :
        continue
    if imp2_orig in aln:
        imp2_tp = aln[imp2_orig]
    else :
        continue
    for imp1 in imp1_tp :
        for imp2 in imp2_tp :
            if (imp1, imp2) in energys :
                print(('Error: (imp1,imp2) in energys = True',imp1_orig,imp2_orig,imp1,imp2))
            energys[(imp1,imp2)] = float(linesp[15])
f_edc.close()
    
for con in ninfo.contacts:
    if (con.imp1, con.imp2) in energys :
        e = energys[(con.imp1, con.imp2)]
        if e > 0.0 :
            con.coef = 0.0
        else :
            con.coef = factor * (-1.0 * e) * con.coef
        
f_out_ninfo.write_all(ninfo)
f_out_ninfo.close()    
