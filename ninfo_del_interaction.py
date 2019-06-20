#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 2011/05/05 coded by Naoto HORI

import sys
from file_ninfo import NinfoFile
from ninfo import NinfoSet

if len(sys.argv) != 4:
    print() 
    print('Usage: % SCRIPT [NINFO file] [command file]  [output NINFO file]')
    print() 
    sys.exit(2)

# Input & file opening
f_ninfo = NinfoFile(sys.argv[1])
f_ninfo.open_to_read()

f_cmd = open(sys.argv[2], 'r')

f_out = NinfoFile(sys.argv[3])
f_out.open_to_write()

# Read ninfo file
ns = NinfoSet()
f_ninfo.read_all(ns)
f_ninfo.close()

# Read command file
cmd_del = []
for line in f_cmd :
    if line[:3] == 'del' :
        cmd_del.append(line.split())
        
# delete interaction
for icmd in cmd_del :
    cmd = icmd[0]
    type = icmd[1]
    target_1 = int(icmd[2])
    if icmd[3] in ('*','!') :
        target_2 = icmd[3]
    else:
        target_2 = int(icmd[3])
    
    if cmd != 'del' :
        continue
    
    if type in ('contact','any') :
        if target_2 == '*' :
            ns.contacts = [
               con for con in ns.contacts 
               if not (con.iunit1==target_1 or con.iunit2==target_1)
               ]
        elif  target_2 == '!' :
            ns.contacts = [
               con for con in ns.contacts 
               if not (    (con.iunit1==target_1 and con.iunit2!=target_1) 
                        or (con.iunit1!=target_1 and con.iunit2==target_1) )
               ]
        else :
            ns.contacts = [
               con for con in ns.contacts 
               if not (    (con.iunit1==target_1 and con.iunit2==target_2) 
                        or (con.iunit1==target_2 and con.iunit2==target_1) )
               ]
            
    if type in ('basepair','any') :
        if target_2 == '*' :
            ns.basepairs = [
               bp for bp in ns.basepairs 
               if not (bp.iunit1==target_1 or bp.iunit2==target_1)
               ]
        elif  target_2 == '!' :
            ns.basepairs = [
               bp for bp in ns.basepairs 
               if not (    (bp.iunit1==target_1 and bp.iunit2!=target_1) 
                        or (bp.iunit1!=target_1 and bp.iunit2==target_1) )
               ]
        else :
            ns.basepairs = [
               bp for bp in ns.basepairs 
               if not (    (bp.iunit1==target_1 and bp.iunit2==target_2) 
                        or (bp.iunit1==target_2 and bp.iunit2==target_1) )
               ]

# Output
ns.reassign_id()
f_out.write_all(ns)
f_out.close()

'''
== command file ==
del bp 1 *
del contact 2 3
del any 4 *
del basepair 5 6
'''