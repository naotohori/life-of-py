#!/usr/bin/env python
'''
@author: Naoto Hori
'''

import sys
from cafysis.file_io.ninfo import NinfoFile
from cafysis.elements.ninfo import NinfoSet

if len(sys.argv) < 4 :
    print "Usage: SCRIPT [input ninfo] [factor] [unit ID ...] [output ninfo]"
    sys.exit(2)
    
f_ninfo_in = NinfoFile(sys.argv[1])
f_ninfo_in.open_to_read()
ninfo = NinfoSet()
f_ninfo_in.read_all(ninfo)
f_ninfo_in.close()

k_factor = float(sys.argv[2])

f_ninfo_out = NinfoFile(sys.argv[-1])
f_ninfo_out.open_to_write()

ids = []
for arg in sys.argv[3:-1] :
    ids.append(int(arg))
#print ids

for con in ninfo.contacts :
    # intra-chain is not treated
    if con.iunit1 == con.iunit2 :
        continue
    
    if con.iunit1 in ids or con.iunit2 in ids:
        con.coef = con.coef * k_factor * con.factor
        con.factor = 1.0
        
f_ninfo_out.write_all(ninfo)
f_ninfo_out.close()        