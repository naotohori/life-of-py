#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@author: Naoto Hori
'''

import sys
from cafysis.file_io.ninfo import NinfoFile, NinfoSet

def show_usage() :
    print ''
    print ' usage: ninfo_split.py [ninfo file] [prefix (can be dir)]'
    print ''
    
if len(sys.argv) != 3:
    show_usage()
    sys.exit(2)
    
f_ninfo = NinfoFile(sys.argv[1])
f_ninfo.open_to_read()
ninfo = NinfoSet()
f_ninfo.read_all(ninfo)
f_ninfo.close()

out_prefix = sys.argv[2]
ninfo.update_info()
nUnit = ninfo.max_unit
print nUnit 
ninfo_out = {}
for i in xrange(1,nUnit+1) :
    for j in xrange(i,nUnit+1) :
        ninfo_out[(i,j)] = NinfoSet()
        
for bd in ninfo.bondlengths :
    ninfo_out[(bd.iunit1, bd.iunit2)].bondlengths.append(bd)
    
for ba in ninfo.bondangles :
    ninfo_out[(ba.iunit1, ba.iunit2)].bondangles.append(ba)
    
for dih in ninfo.dihedrals :
    ninfo_out[(dih.iunit1, dih.iunit2)].dihedrals.append(dih)
    
for con in ninfo.contacts :
    ninfo_out[(con.iunit1, con.iunit2)].contacts.append(con)
    
for bp in ninfo.basepairs:
    ninfo_out[(bp.iunit1, bp.iunit2)].basepairs.append(bp)
    
for bs in ninfo.basestacks :
    ninfo_out[(bs.iunit1, bs.iunit2)].basestacks.append(bs)
    
for i in xrange(1,nUnit+1) :
    for j in xrange(i, nUnit+1) :
        if (len(ninfo_out[(i,j)].bondlengths) != 0 or len(ninfo_out[(i,j)].bondangles) != 0  or
            len(ninfo_out[(i,j)].dihedrals) != 0   or len(ninfo_out[(i,j)].contacts) != 0    or
            len(ninfo_out[(i,j)].basepairs) != 0    or len(ninfo_out[(i,j)].basestacks) != 0 ):
            filename = out_prefix + '%i_%i' % (i,j) + '.ninfo'
            nf = NinfoFile(filename)
            nf.open_to_write()
            nf.write_all(ninfo_out[(i,j)])