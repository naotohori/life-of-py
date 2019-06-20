#!/usr/bin/env python

from file_ninfo import NinfoFile
from ninfo import NinfoSet
from copy import deepcopy
import sys

def show_usage() :
    print() 
    print('Usage: % SCRIPT [ninfo (all)] [ninfo (part)] [ninfo (output)] [[(unit ID) (starting imp)] ....]')
    print() 
    
if len(sys.argv) < 5:
    show_usage()
    sys.exit(2)
    
# input
f_ninfo_base = NinfoFile(sys.argv[1])
f_ninfo_part = NinfoFile(sys.argv[2])
f_ninfo_out  = NinfoFile(sys.argv[3])
target_units = []
target_imp_start = []

n = 1
for arg in sys.argv[4:] :
    if n == 1:
        target_units.append(int(arg))
    else :
        target_imp_start.append(int(arg))
    n *= -1
if n == -1:
    show_usage()
    sys.exit(2)
    
# for test
#f_ninfo_base = NinfoFile("../test/ribo_005.ninfo")
#f_ninfo_part = NinfoFile("../test/3l0u_fmat012.ninfo")
#f_ninfo_out = NinfoFile("../test/new.ninfo")
#target_units = [52,53,54]
#target_imp_start = [18902,19129,19356]

# read ninfo(all)
f_ninfo_base.open_to_read()
ninfo_base = NinfoSet()
f_ninfo_base.read_all(ninfo_base)
f_ninfo_base.close()

# read ninfo(part)
f_ninfo_part.open_to_read()
ninfo_part = NinfoSet()
f_ninfo_part.read_all(ninfo_part)
f_ninfo_part.close()

# delete
for target_id in target_units :
    ninfo_base.bondlengths = [bd for bd in ninfo_base.bondlengths if not (bd.iunit1==target_id and bd.iunit2==target_id)]
    ninfo_base.bondangles = [ba for ba in ninfo_base.bondangles if not (ba.iunit1==target_id and ba.iunit2==target_id)]
    ninfo_base.dihedrals = [dih for dih in ninfo_base.dihedrals if not (dih.iunit1==target_id and dih.iunit2==target_id)]
    ninfo_base.contacts = [con for con in ninfo_base.contacts if not (con.iunit1==target_id and con.iunit2==target_id)]
    ninfo_base.basepairs = [bp for bp in ninfo_base.basepairs if not (bp.iunit1==target_id and bp.iunit2==target_id)]
    ninfo_base.basestacks = [bs for bs in ninfo_base.basestacks if not (bs.iunit1==target_id and bs.iunit2==target_id)]
    
# add
for (iun, target_id) in enumerate(target_units) :
    
    for bd_part in ninfo_part.bondlengths :
        bd_append = deepcopy(bd_part)
        bd_append.iunit1 = target_id
        bd_append.iunit2 = target_id
        bd_append.imp1 = bd_part.imp1 + target_imp_start[iun] - 1
        bd_append.imp2 = bd_part.imp2 + target_imp_start[iun] - 1
        ninfo_base.bondlengths.append(bd_append)
        
    for ba_part in ninfo_part.bondangles :
        ba_append = deepcopy(ba_part)
        ba_append.iunit1 = target_id
        ba_append.iunit2 = target_id
        ba_append.iunit3 = target_id
        ba_append.imp1 = ba_part.imp1 + target_imp_start[iun] - 1
        ba_append.imp2 = ba_part.imp2 + target_imp_start[iun] - 1
        ba_append.imp3 = ba_part.imp3 + target_imp_start[iun] - 1
        ninfo_base.bondangles.append(ba_append)
        
    for dih_part in ninfo_part.dihedrals :
        dih_append = deepcopy(dih_part)
        dih_append.iunit1 = target_id
        dih_append.iunit2 = target_id
        dih_append.iunit3 = target_id
        dih_append.iunit4 = target_id
        dih_append.imp1 = dih_part.imp1 + target_imp_start[iun] - 1
        dih_append.imp2 = dih_part.imp2 + target_imp_start[iun] - 1
        dih_append.imp3 = dih_part.imp3 + target_imp_start[iun] - 1
        dih_append.imp4 = dih_part.imp4 + target_imp_start[iun] - 1
        ninfo_base.dihedrals.append(dih_append)
        
    for con_part in ninfo_part.contacts :
        con_append = deepcopy(con_part)
        con_append.iunit1 = target_id
        con_append.iunit2 = target_id
        con_append.imp1 = con_part.imp1 + target_imp_start[iun] - 1
        con_append.imp2 = con_part.imp2 + target_imp_start[iun] - 1
        ninfo_base.contacts.append(con_append)
        
    for bp_part in ninfo_part.basepairs :
        bp_append = deepcopy(bp_part)
        bp_append.iunit1 = target_id
        bp_append.iunit2 = target_id
        bp_append.imp1 = bp_part.imp1 + target_imp_start[iun] - 1
        bp_append.imp2 = bp_part.imp2 + target_imp_start[iun] - 1
        ninfo_base.basepairs.append(bp_append)
        
    for bs_part in ninfo_part.basestacks :
        bs_append = deepcopy(bs_part)
        bs_append.iunit1 = target_id
        bs_append.iunit2 = target_id
        bs_append.imp1 = bs_part.imp1 + target_imp_start[iun] - 1
        bs_append.imp2 = bs_part.imp2 + target_imp_start[iun] - 1
        ninfo_base.basestacks.append(bs_append)
        
# re-assign id
ninfo_base.reassign_id()

f_ninfo_out.open_to_write()
f_ninfo_out.write_all(ninfo_base)
f_ninfo_out.close()
