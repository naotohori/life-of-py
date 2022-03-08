#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2011/12/14
@author: Naoto Hori
'''

if __name__ == '__main__':
    import sys
    if not len(sys.argv) in (4,) :
        print('Usage: % SCRIPT [ninfo file(i)] [cafemol PDB file(i)] [ptraj file(o)]')
        sys.exit(2)
    
    path_output = 'measure/'
    time = 1.0
    digit_for_id = 3
    def_pro = "CA"
    def_P = "P"
    def_S = "C4',O4',C1',C2',C3'"
    def_RAG = "N1"
    def_YUC = "N3"
    
    
    from lop.file_io.ninfo import NinfoFile
    from lop.elements.ninfo import NinfoSet
    file_ninfo = NinfoFile(sys.argv[1]) 
    file_ninfo.open_to_read()
    ns = NinfoSet()
    file_ninfo.read_all(ns)
    
    from lop.file_io.pdb import PdbFile
    file_pdb = PdbFile(sys.argv[2])
    file_pdb.open_to_read()
    chains = file_pdb.read_all()
    
    imp2ires = []
    ires_total = 0
    ires_prev = -1
    for c in chains:
        for i in range(c.num_atom()):
            ires = c.get_atom(i).res_seq
            if ires != ires_prev:
                ires_total += 1
                ires_prev = ires
            imp2ires.append(ires_total)
    
    file_ptraj = open(sys.argv[-1],'w')
    print(imp2ires)
    print(len(imp2ires))
    
    for bd in ns.bondlengths:
        type_out = bd.type
        if bd.type == 'pp':
            com1 = def_pro
            com2 = def_pro
        elif bd.type == 'PS':
            com1 = def_P
            com2 = def_S
        elif bd.type in ('SR','SA','SG'):
            com1 = def_S
            com2 = def_RAG
            type_out = 'SR'
        elif bd.type in ('SY','SU','SC'):
            com1 = def_S
            com2 = def_YUC
            type_out = 'SY'
        elif bd.type == 'SP':
            com1 = def_S
            com2 = def_P
        else:
            print('Error: invalid bond type, '+bd.type)
            sys.exit(2)
            
        cmd = ('distance bond_%i_%i :%i@%s :%i@%s out %s/bond_%s_%0'
             + ('%i' % digit_for_id) + 'i_%0' + ('%i' % digit_for_id)
             + 'i.out time %f\n')
        ires1 = imp2ires[bd.imp1-1]
        ires2 = imp2ires[bd.imp2-1]
        file_ptraj.write(cmd % (bd.imp1, bd.imp2, ires1, com1, ires2, com2,
                                path_output, type_out, bd.imp1, bd.imp2, time) )
        
    for ba in ns.bondangles:
        type_out = ba.type
        if ba.type == 'ppp':
            com1 = def_pro
            com2 = def_pro
            com3 = def_pro
        elif ba.type in ('PSR','PSA','PSG'):
            com1 = def_P
            com2 = def_S
            com3 = def_RAG
            type_out = 'PSR'
        elif ba.type in ('PSY','PSU','PSC'):
            com1 = def_P
            com2 = def_S
            com3 = def_YUC
            type_out = 'PSY'
        elif ba.type in ('RSP','ASP','GSP'):
            com1 = def_RAG
            com2 = def_S
            com3 = def_P
            type_out = 'RSP'
        elif ba.type in ('YSP','USP','CSP'):
            com1 = def_YUC
            com2 = def_S
            com3 = def_P
            type_out = 'YSP'
        elif ba.type == 'PSP':
            com1 = def_P
            com2 = def_S
            com3 = def_P
        elif ba.type == 'SPS':
            com1 = def_S
            com2 = def_P
            com3 = def_S
        else:
            print('Error: invalid angle type, '+ba.type)
            sys.exit(2)
        
        cmd = ('angle angl_%i_%i_%i :%i@%s :%i@%s :%i@%s out %s/angl_%s_%0'
             + ('%i' % digit_for_id) + 'i_%0' + ('%i' % digit_for_id) + 'i_%0'
             + ('%i' % digit_for_id) + 'i.out time %f\n')
        ires1 = imp2ires[ba.imp1-1]
        ires2 = imp2ires[ba.imp2-1]
        ires3 = imp2ires[ba.imp3-1]
        file_ptraj.write(cmd % (ba.imp1, ba.imp2, ba.imp3,
                                ires1, com1, ires2, com2, ires3, com3,
                                path_output, type_out,
                                ba.imp1, ba.imp2, ba.imp3, time) )
    
    for dih in ns.dihedrals:
        type_out = dih.type
        if dih.type == 'pppp':
            com1 = def_pro
            com2 = def_pro
            com3 = def_pro
            com4 = def_pro
        elif dih.type in ('SPSR','SPSA','SPSG'):
            com1 = def_S
            com2 = def_P
            com3 = def_S
            com4 = def_RAG
            type_out = 'SPSR'
        elif dih.type in ('SPSY','SPSU','SPSC'):
            com1 = def_S
            com2 = def_P
            com3 = def_S
            com4 = def_YUC
            type_out = 'SPSY'
        elif dih.type in ('RSPS','ASPS','GSPS'):
            com1 = def_RAG
            com2 = def_S
            com3 = def_P
            com4 = def_S
            type_out = 'RSPS'
        elif dih.type in ('YSPS','USPS','CSPS'):
            com1 = def_YUC
            com2 = def_S
            com3 = def_P
            com4 = def_S
            type_out = 'YSPS'
        elif dih.type == 'SPSP':
            com1 = def_S
            com2 = def_P
            com3 = def_S
            com4 = def_P
        elif dih.type == 'PSPS':
            com1 = def_P
            com2 = def_S
            com3 = def_P
            com4 = def_S
        elif dih.type[1:3] == 'SS':
            continue
        else:
            print('Error: invalid dihedral type, ' + dih.type)
            sys.exit(2)
            
        cmd = ('dihedral dih_%i_%i_%i_%i :%i@%s :%i@%s :%i@%s :%i@%s out %s/dih_%s_%0'
             + ('%i' % digit_for_id) + 'i_%0' + ('%i' % digit_for_id) + 'i_%0'
             + ('%i' % digit_for_id) + 'i_%0' + ('%i' % digit_for_id) 
             + 'i.out time %f\n')
        ires1 = imp2ires[dih.imp1-1]
        ires2 = imp2ires[dih.imp2-1]
        ires3 = imp2ires[dih.imp3-1]
        ires4 = imp2ires[dih.imp4-1]
        file_ptraj.write(cmd % (dih.imp1, dih.imp2, dih.imp3, dih.imp4,
                                ires1, com1, ires2, com2, ires3, com3, ires4, com4,
                                path_output, type_out, 
                                dih.imp1, dih.imp2, dih.imp3, dih.imp4, time) )
        
    for con in ns.contacts:
        type_out = con.type
        if con.type[0] == 'p':
            com1 = def_pro
        elif con.type[0] == 'P':
            com1 = def_P
        elif con.type[0] == 'S':
            com1 = def_S
        elif con.type[0] in ('A','G','R'):
            com1 = def_RAG
        elif con.type[0] in ('Y','U','C'):
            com1 = def_YUC
        else:
            print('Error: invalid contact type, '+con.type)
            sys.exit(2)
        if con.type[2] == 'p':
            com2 = def_pro
        elif con.type[2] == 'P':
            com2 = def_P
        elif con.type[2] == 'S':
            com2 = def_S
        elif con.type[2] in ('A','G','R'):
            com2 = def_RAG
        elif con.type[2] in ('Y','U','C'):
            com2 = def_YUC
        else:
            print('Error: invalid contact type, '+con.type)
            sys.exit(2)
        if type_out in ('p-A','p-G','p-R','p-Y','p-U','p-C'):
            type_out = 'p-B'
        elif type_out in ('P-A','P-G','P-R','P-Y','P-U','P-C'):
            type_out = 'P-B'
        elif type_out in ('S-A','S-G','S-R','S-Y','S-U','S-C'):
            type_out = 'S-B'
        elif (type_out[0] in ('A','G','R','Y','U','C') and 
              type_out[2] in ('A','G','R','Y','U','C')):
            type_out = 'B-B'
            
        cmd = ('distance con_%i_%i :%i@%s :%i@%s out %s/con_%s_%0'
             + ('%i' % digit_for_id) + 'i_%0' + ('%i' % digit_for_id) 
             + 'i.out time %f\n')
        ires1 = imp2ires[con.imp1-1]
        ires2 = imp2ires[con.imp2-1]
        file_ptraj.write(cmd % (con.imp1, con.imp2,
                                ires1, com1, ires2, com2,
                                path_output, type_out,
                                con.imp1, con.imp2, time) )
    
    for bp in ns.basepairs:
        if bp.type[0] in ('A','G','R'):
            com1 = def_RAG
        elif bp.type[0] in ('Y','U','C'):
            com1 = def_YUC
        else:
            print('Error: invalid basepair type, '+bp.type)
            sys.exit(2)
        if bp.type[2] in ('A','G','R'):
            com2 = def_RAG
        elif bp.type[2] in ('Y','U','C'):
            com2 = def_YUC
        else:
            print('Error: invalid basepair type, '+bp.type)
            sys.exit(2)
            
        cmd = ('distance bp_%i_%i :%i@%s :%i@%s out %s/bp_HB%i_%0'
             + ('%i' % digit_for_id) + 'i_%0' + ('%i' % digit_for_id) 
             + 'i.out time %f\n')
        ires1 = imp2ires[bp.imp1-1]
        ires2 = imp2ires[bp.imp2-1]
        file_ptraj.write(cmd % (bp.imp1, bp.imp2,
                                ires1, com1, ires2, com2,
                                path_output, bp.nhb, bp.imp1, bp.imp2, time) )
    
    for bs in ns.basestacks:
        if bs.type[0] in ('A','G','R'):
            com1 = def_RAG
        elif bs.type[0] in ('Y','U','C'):
            com1 = def_YUC
        else:
            print('Error: invalid basestack type, '+bs.type)
            sys.exit(2)
        if bs.type[2] in ('A','G','R'):
            com2 = def_RAG
        elif bs.type[2] in ('Y','U','C'):
            com2 = def_YUC
        else:
            print('Error: invalid basestack type, '+bs.type)
            sys.exit(2)
            
        cmd = ('distance bs_%i_%i :%i@%s :%i@%s out %s/bs_%0'
             + ('%i' % digit_for_id) + 'i_%0' + ('%i' % digit_for_id) 
             + 'i.out time %f\n')
        ires1 = imp2ires[bs.imp1-1]
        ires2 = imp2ires[bs.imp2-1]
        file_ptraj.write(cmd % (bs.imp1, bs.imp2,
                                ires1, com1, ires2, com2,
                                path_output, bs.imp1, bs.imp2, time) )
