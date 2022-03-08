#!/usr/bin/env python

import copy
from lop.file_io.ninfo import NinfoFile
from lop.elements.ninfo import NinfoSet

Nunit_per_repeat = 2
Nmp_per_repeat = 130   # The last mp of the unit dimer (dsRNA)
Nrepeat = 27    # Number of repeat 

ns = NinfoSet()

ns_file = NinfoFile('./miR21inh_noNN.ninfo')
ns_file.open_to_read()
ns_file.read_all(ns)
ns_file.close()

ns_new = NinfoSet()

for irepeat in range(Nrepeat):

    unit_offset = irepeat * Nunit_per_repeat
    mp_offset = irepeat * Nmp_per_repeat

    for ni in ns.bondlengths:
        ni_new = copy.deepcopy(ni)

        ni_new.iunit1 += unit_offset
        ni_new.iunit2 += unit_offset
        ni_new.imp1 += mp_offset
        ni_new.imp2 += mp_offset
        
        ns_new.bondlengths.append(ni_new)

    for ni in ns.bondangles:
        ni_new = copy.deepcopy(ni)

        ni_new.iunit1 += unit_offset
        ni_new.iunit2 += unit_offset
        ni_new.imp1 += mp_offset
        ni_new.imp2 += mp_offset
        ni_new.imp3 += mp_offset
        
        ns_new.bondangles.append(ni_new)

    for ni in ns.basestackDTs:
        ni_new = copy.deepcopy(ni)

        ni_new.iunit1 += unit_offset
        ni_new.iunit2 += unit_offset
        ni_new.imp1 += mp_offset
        ni_new.imp2 += mp_offset
        ni_new.dih1_iunit1 += unit_offset
        ni_new.dih1_iunit2 += unit_offset
        ni_new.dih1_imp1 += mp_offset
        ni_new.dih1_imp2 += mp_offset
        ni_new.dih1_imp3 += mp_offset
        ni_new.dih1_imp4 += mp_offset
        ni_new.dih2_iunit1 += unit_offset
        ni_new.dih2_iunit2 += unit_offset
        ni_new.dih2_imp1 += mp_offset
        ni_new.dih2_imp2 += mp_offset
        ni_new.dih2_imp3 += mp_offset
        ni_new.dih2_imp4 += mp_offset
        
        ns_new.basestackDTs.append(ni_new)

    for ni in ns.hbondDTs:
        ni_new = copy.deepcopy(ni)

        ni_new.iunit1 += unit_offset
        ni_new.iunit2 += unit_offset
        ni_new.imp1 += mp_offset
        ni_new.imp2 += mp_offset
        ni_new.ang1_iunit1 += unit_offset
        ni_new.ang1_iunit2 += unit_offset
        ni_new.ang1_imp1 += mp_offset
        ni_new.ang1_imp2 += mp_offset
        ni_new.ang1_imp3 += mp_offset
        ni_new.ang2_iunit1 += unit_offset
        ni_new.ang2_iunit2 += unit_offset
        ni_new.ang2_imp1 += mp_offset
        ni_new.ang2_imp2 += mp_offset
        ni_new.ang2_imp3 += mp_offset
        ni_new.dih0_iunit1 += unit_offset
        ni_new.dih0_iunit2 += unit_offset
        ni_new.dih0_imp1 += mp_offset
        ni_new.dih0_imp2 += mp_offset
        ni_new.dih0_imp3 += mp_offset
        ni_new.dih0_imp4 += mp_offset
        ni_new.dih1_iunit1 += unit_offset
        ni_new.dih1_iunit2 += unit_offset
        ni_new.dih1_imp1 += mp_offset
        ni_new.dih1_imp2 += mp_offset
        ni_new.dih1_imp3 += mp_offset
        ni_new.dih1_imp4 += mp_offset
        ni_new.dih2_iunit1 += unit_offset
        ni_new.dih2_iunit2 += unit_offset
        ni_new.dih2_imp1 += mp_offset
        ni_new.dih2_imp2 += mp_offset
        ni_new.dih2_imp3 += mp_offset
        ni_new.dih2_imp4 += mp_offset
        
        ns_new.hbondDTs.append(ni_new)

    for ni in ns.tertiarystackDTs:
        ni_new = copy.deepcopy(ni)

        ni_new.iunit1 += unit_offset
        ni_new.iunit2 += unit_offset
        ni_new.imp1 += mp_offset
        ni_new.imp2 += mp_offset
        ni_new.ang1_iunit1 += unit_offset
        ni_new.ang1_iunit2 += unit_offset
        ni_new.ang1_imp1 += mp_offset
        ni_new.ang1_imp2 += mp_offset
        ni_new.ang1_imp3 += mp_offset
        ni_new.ang2_iunit1 += unit_offset
        ni_new.ang2_iunit2 += unit_offset
        ni_new.ang2_imp1 += mp_offset
        ni_new.ang2_imp2 += mp_offset
        ni_new.ang2_imp3 += mp_offset
        ni_new.dih0_iunit1 += unit_offset
        ni_new.dih0_iunit2 += unit_offset
        ni_new.dih0_imp1 += mp_offset
        ni_new.dih0_imp2 += mp_offset
        ni_new.dih0_imp3 += mp_offset
        ni_new.dih0_imp4 += mp_offset
        ni_new.dih1_iunit1 += unit_offset
        ni_new.dih1_iunit2 += unit_offset
        ni_new.dih1_imp1 += mp_offset
        ni_new.dih1_imp2 += mp_offset
        ni_new.dih1_imp3 += mp_offset
        ni_new.dih1_imp4 += mp_offset
        ni_new.dih2_iunit1 += unit_offset
        ni_new.dih2_iunit2 += unit_offset
        ni_new.dih2_imp1 += mp_offset
        ni_new.dih2_imp2 += mp_offset
        ni_new.dih2_imp3 += mp_offset
        ni_new.dih2_imp4 += mp_offset
        
        ns_new.tertiarystackDTs.append(ni_new)

ns_new.update_info()
ns_new.reassign_id()


ns_file = NinfoFile('./new.ninfo')
ns_file.open_to_write()
ns_file.write_all(ns_new)
ns_file.close()
