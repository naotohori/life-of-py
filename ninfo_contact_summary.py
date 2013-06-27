#!/usr/bin/env python
'''
@author: Naoto Hori
'''

from cafysis.file_io.ninfo import NinfoFile
from cafysis.elements.ninfo import NinfoSet
import sys

if not len(sys.argv) in (3,4):
    print ('\nUsage 1: % SCRIPT [ninfo file] [output file (csv)]')
    print ('      2: % SCRIPT [ninfo file] [output file (csv)] [output list file]\n')
    sys.exit(2)
    
# input
f_ninfo = NinfoFile(sys.argv[1])
f_out = open(sys.argv[2],'w')
flg_out_list = False
if len(sys.argv) == 4 :
    flg_out_list = True
    f_out_list = open(sys.argv[3],'w')

ninfo = NinfoSet()
f_ninfo.open_to_read()
f_ninfo.read_all(ninfo)
f_ninfo.close()

ninfo.update_info()
nunit = ninfo.max_unit

# Count the number of contacts or basepairs
data = {} 
for i in xrange(1,nunit+1) :
    for j in xrange(i, nunit+1) :
        n_con = len(ninfo.get_contacts_by_unit(i, j))
        n_bp  = len(ninfo.get_basepairs_by_unit(i, j))
        data[(i,j)] = (n_con, n_bp)
        data[(j,i)] = (n_con, n_bp)

####### Output for CSV
### contact
# First line
f_out.write("#contact")
for j in xrange(1, nunit+1) :
    f_out.write(",%i" % j)
f_out.write("\n")

# Data
for i in xrange(1, nunit+1) :
    f_out.write("%i" % i)
    for j in xrange(1, nunit+1) :
        f_out.write(",%i" % data[(i,j)][0] )
    f_out.write("\n")
    
### basepair
# First line
f_out.write("\n\n#basepair")
for j in xrange(1, nunit+1) :
    f_out.write(",%i" % j)
f_out.write("\n")

# Data
for i in xrange(1, nunit+1) :
    f_out.write("%i" % i)
    for j in xrange(1, nunit+1) :
        f_out.write(",%i" % data[(i,j)][1] )
    f_out.write("\n")
f_out.close()    

####### Output for list
if flg_out_list :
    for i in xrange(1, nunit+1) :
        for j in xrange(i, nunit+1) :
            f_out_list.write("%5i %5i %10i %10i\n" % ((i,j) + data[(i,j)]) )
        f_out_list.write('#--------------------------------\n')
    f_out_list.close()
    
