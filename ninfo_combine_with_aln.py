#!/usr/bin/env python

def show_usage() :
    print('')
    print(' usage: SCRIPT [ninfo file 1] [ninfo file 2] [cmd file] [ninfo output file]')
    print('')
    
import sys
from file_ninfo import NinfoFile
from ninfo import NinfoSet
from copy import deepcopy

if len(sys.argv) != 5 :
    show_usage()
    sys.exit(2)
    
ninfo1 = NinfoSet()
f_ninfo = NinfoFile(sys.argv[1])
f_ninfo.open_to_read()
f_ninfo.read_all(ninfo1)
f_ninfo.close()

ninfo2 = NinfoSet()
f_ninfo = NinfoFile(sys.argv[2])
f_ninfo.open_to_read()
f_ninfo.read_all(ninfo2)
f_ninfo.close()

f_out = NinfoFile(sys.argv[-1])
f_out.open_to_write()

f_cmd = open(sys.argv[3], 'r')
cmds = []
for line in f_cmd :
    if line.startswith(' ') or line.startswith('#') or line.startswith('*') :
        continue
    linesp = line.split()
    # 40_55.ninfo 40 43 L25.ninfo 55 54 AtRNA.ninfo 
    cmds.append( (linesp[0],   # ninfo file 2
                  int(linesp[1]), int(linesp[2]), linesp[3],   # chainID in 2, chainID in 1, aln
                  int(linesp[4]), int(linesp[5]), linesp[6]) ) # chainID in 2, chainID in 1, aln
    

for cmd in cmds :
    if cmd[0] == 'input' :
        subset = ninfo2
    else :
        subset = NinfoSet()
        f = NinfoFile(cmd[0])
        f.open_to_read()
        f.read_all(subset)
        f.close()
    sub_iunit1 = cmd[1]
    main_iunit1 = cmd[2]
    filename_aln1 = cmd[3]
    sub_iunit2 = cmd[4]
    main_iunit2 = cmd[5]
    filename_aln2 = cmd[6]
    
    # aln1
    aln1 = {}
    f_aln = open(filename_aln1, 'r')
    for line in f_aln :
        if line.find('#') != -1 :
            continue
        linesp = line.split()
        imp1 = int(linesp[1])
        imp2 = int(linesp[3])
        #dist = float(linesp[4])
        aln1[imp1] = imp2
    f_aln.close()
    
    # aln2
    aln2 = {}
    f_aln = open(filename_aln2, 'r')
    for line in f_aln :
        if line.find('#') != -1 :
            continue
        linesp = line.split()
        imp1 = int(linesp[1])
        imp2 = int(linesp[3])
        #dist = float(linesp[4])
        aln2[imp1] = imp2
    f_aln.close()
        
    for con in subset.contacts :
        if con.iunit1 != sub_iunit1 or con.iunit2 != sub_iunit2 :
            continue
        con_append = deepcopy(con)
        if con.imp1 in aln1 :
            con_append.imp1 = aln1[con.imp1]
        else :
            continue
        if con.imp2 in aln2 :
            con_append.imp2 = aln2[con.imp2]
        else :
            continue
        con_append.iunit1 = main_iunit1
        con_append.iunit2 = main_iunit2
        ninfo1.contacts.append(con_append)
        
    for bp in subset.basepairs :
        if bp.iunit1 != sub_iunit1 or bp.iunit2 != sub_iunit2 :
            continue
        bp_append = deepcopy(bp)
        if bp.imp1 in aln1 :
            bp_append.imp1 = aln1[bp.imp1]
        else :
            continue
        if bp.imp2 in aln2 :
            bp_append.imp2 = aln2[bp.imp2]
        else :
            continue
        bp_append.iunit1 = main_iunit1
        bp_append.iunit2 = main_iunit2
        ninfo1.basepairs.append(bp_append)
        
ninfo1.reassign_id()

f_out.write_all(ninfo1)
f_out.close()