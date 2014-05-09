#!/usr/bin/env python

import sys
from file_ninfo import NinfoFile
from ninfo import NinfoSet, Contact, BasePair

if len(sys.argv) != 4:
    print ('Usage: SCRIPT [ninfo file] [tRNA information file] [output ninfo]')
    sys.exit(2)
    
f_ninfo = NinfoFile(sys.argv[1])
ninfo = NinfoSet()
f_ninfo.open_to_read()
f_ninfo.read_all(ninfo)
f_ninfo.close()

f_trna = open(sys.argv[2], 'r')
for line in f_trna :
    if line.find('#') != -1 :
        continue
    linesp = line.split()
    if linesp[0] == 'A' :
        unit_A = int(linesp[1])
        imp_pair_A = (int(linesp[2]), int(linesp[3]))
    elif linesp[0] == 'P' :
        unit_P = int(linesp[1])
        imp_pair_P = (int(linesp[2]), int(linesp[3]))
    elif linesp[0] == 'E' :
        unit_E = int(linesp[1])
        imp_pair_E = (int(linesp[2]), int(linesp[3]))
    elif linesp[0] == 'm' :
        unit_m = int(linesp[1])
        imp_pair_m = (int(linesp[2]), int(linesp[3]))
        
if unit_A == None or unit_P == None or unit_E == None or unit_m == None :
    print 'Error:  tRNA information is not enough'
    sys.exit(2)
        
f_out_ninfo = NinfoFile(sys.argv[3])
f_out_ninfo.open_to_write()
    
new_contacts = []
for con in ninfo.contacts :
    if con.iunit1 == con.iunit2 :
        continue
    if (  (con.iunit1 == unit_A and con.iunit2 != unit_P and con.iunit2 != unit_E and con.iunit2 != unit_m) 
        or(con.iunit2 == unit_A and con.iunit1 != unit_P and con.iunit1 != unit_E and con.iunit1 != unit_m)) :
        newconP = Contact()
        newconE = Contact()
        if con.iunit1 == unit_A :
            newconP.imp1 = imp_pair_P[0] + (con.imp1 - imp_pair_A[0])
            newconE.imp1 = imp_pair_E[0] + (con.imp1 - imp_pair_A[0])
            newconP.imp2 = con.imp2
            newconE.imp2 = con.imp2
            newconP.iunit1 = unit_P
            newconE.iunit1 = unit_E
            newconP.iunit2 = con.iunit2
            newconE.iunit2 = con.iunit2
        else :
            newconP.imp2 = imp_pair_P[0] + (con.imp2 - imp_pair_A[0])
            newconE.imp2 = imp_pair_E[0] + (con.imp2 - imp_pair_A[0])
            newconP.imp1 = con.imp1
            newconE.imp1 = con.imp1
            newconP.iunit2 = unit_P
            newconE.iunit2 = unit_E
            newconP.iunit1 = con.iunit1
            newconE.iunit1 = con.iunit1
        newconP.imp1un = con.imp1un
        newconE.imp1un = con.imp1un
        newconP.imp2un = con.imp2un
        newconE.imp2un = con.imp2un
        newconP.type = con.type
        newconE.type = con.type
        newconP.coef = con.coef
        newconE.coef = con.coef
        newconP.factor = con.factor
        newconE.factor = con.factor
        newconP.dummy = con.dummy
        newconE.dummy = con.dummy
        newconP.correct_mgo = con.correct_mgo
        newconE.correct_mgo = con.correct_mgo
        newconP.native = con.native
        newconE.native = con.native
        new_contacts.append(newconP)
        new_contacts.append(newconE)
        
    if (  (con.iunit1 == unit_P and con.iunit2 != unit_A and con.iunit2 != unit_E and con.iunit2 != unit_m) 
        or(con.iunit2 == unit_P and con.iunit1 != unit_A and con.iunit1 != unit_E and con.iunit1 != unit_m)) :
        newconA = Contact()
        newconE = Contact()
        if con.iunit1 == unit_P :
            newconA.imp1 = imp_pair_A[0] + (con.imp1 - imp_pair_P[0])
            newconE.imp1 = imp_pair_E[0] + (con.imp1 - imp_pair_P[0])
            newconA.imp2 = con.imp2
            newconE.imp2 = con.imp2
            newconA.iunit1 = unit_A
            newconE.iunit1 = unit_E
            newconA.iunit2 = con.iunit2
            newconE.iunit2 = con.iunit2
        else :
            newconA.imp2 = imp_pair_A[0] + (con.imp2 - imp_pair_P[0])
            newconE.imp2 = imp_pair_E[0] + (con.imp2 - imp_pair_P[0])
            newconA.imp1 = con.imp1
            newconE.imp1 = con.imp1
            newconA.iunit2 = unit_A
            newconE.iunit2 = unit_E
            newconA.iunit1 = con.iunit1
            newconE.iunit1 = con.iunit1
        newconA.imp1un = con.imp1un
        newconE.imp1un = con.imp1un
        newconA.imp2un = con.imp2un
        newconE.imp2un = con.imp2un
        newconA.type = con.type
        newconE.type = con.type
        newconA.coef = con.coef
        newconE.coef = con.coef
        newconA.factor = con.factor
        newconE.factor = con.factor
        newconA.dummy = con.dummy
        newconE.dummy = con.dummy
        newconA.correct_mgo = con.correct_mgo
        newconE.correct_mgo = con.correct_mgo
        newconA.native = con.native
        newconE.native = con.native
        new_contacts.append(newconA)
        new_contacts.append(newconE)
        
    if (  (con.iunit1 == unit_E and con.iunit2 != unit_A and con.iunit2 != unit_P and con.iunit2 != unit_m) 
        or(con.iunit2 == unit_E and con.iunit1 != unit_A and con.iunit1 != unit_P and con.iunit1 != unit_m)) :
        newconA = Contact()
        newconP = Contact()
        if con.iunit1 == unit_E :
            newconA.imp1 = imp_pair_A[0] + (con.imp1 - imp_pair_E[0])
            newconP.imp1 = imp_pair_P[0] + (con.imp1 - imp_pair_E[0])
            newconA.imp2 = con.imp2
            newconP.imp2 = con.imp2
            newconA.iunit1 = unit_A
            newconP.iunit1 = unit_P
            newconA.iunit2 = con.iunit2
            newconP.iunit2 = con.iunit2
        else :
            newconA.imp2 = imp_pair_A[0] + (con.imp2 - imp_pair_E[0])
            newconP.imp2 = imp_pair_P[0] + (con.imp2 - imp_pair_E[0])
            newconA.imp1 = con.imp1
            newconP.imp1 = con.imp1
            newconA.iunit2 = unit_A
            newconP.iunit2 = unit_P
            newconA.iunit1 = con.iunit1
            newconP.iunit1 = con.iunit1
        newconA.imp1un = con.imp1un
        newconP.imp1un = con.imp1un
        newconA.imp2un = con.imp2un
        newconP.imp2un = con.imp2un
        newconA.type = con.type
        newconP.type = con.type
        newconA.coef = con.coef
        newconP.coef = con.coef
        newconA.factor = con.factor
        newconP.factor = con.factor
        newconA.dummy = con.dummy
        newconP.dummy = con.dummy
        newconA.correct_mgo = con.correct_mgo
        newconP.correct_mgo = con.correct_mgo
        newconA.native = con.native
        newconP.native = con.native
        new_contacts.append(newconA)
        new_contacts.append(newconP)

new_basepairs = []
for con in ninfo.basepairs :
    if con.iunit1 == con.iunit2 :
        continue
    if (  (con.iunit1 == unit_A and con.iunit2 != unit_P and con.iunit2 != unit_E and con.iunit2 != unit_m) 
        or(con.iunit2 == unit_A and con.iunit1 != unit_P and con.iunit1 != unit_E and con.iunit1 != unit_m)) :
        newconP = BasePair()
        newconE = BasePair()
        if con.iunit1 == unit_A :
            newconP.imp1 = imp_pair_P[0] + (con.imp1 - imp_pair_A[0])
            newconE.imp1 = imp_pair_E[0] + (con.imp1 - imp_pair_A[0])
            newconP.imp2 = con.imp2
            newconE.imp2 = con.imp2
            newconP.iunit1 = unit_P
            newconE.iunit1 = unit_E
            newconP.iunit2 = con.iunit2
            newconE.iunit2 = con.iunit2
        else :
            newconP.imp2 = imp_pair_P[0] + (con.imp2 - imp_pair_A[0])
            newconE.imp2 = imp_pair_E[0] + (con.imp2 - imp_pair_A[0])
            newconP.imp1 = con.imp1
            newconE.imp1 = con.imp1
            newconP.iunit2 = unit_P
            newconE.iunit2 = unit_E
            newconP.iunit1 = con.iunit1
            newconE.iunit1 = con.iunit1
        newconP.imp1un = con.imp1un
        newconE.imp1un = con.imp1un
        newconP.imp2un = con.imp2un
        newconE.imp2un = con.imp2un
        newconP.type = con.type
        newconE.type = con.type
        newconP.coef = con.coef
        newconE.coef = con.coef
        newconP.factor = con.factor
        newconE.factor = con.factor
        newconP.dummy = con.dummy
        newconE.dummy = con.dummy
        newconP.correct_mgo = con.correct_mgo
        newconE.correct_mgo = con.correct_mgo
        newconP.native = con.native
        newconE.native = con.native
        newconP.nhb = con.nhb
        newconE.nhb = con.nhb
        new_basepairs.append(newconP)
        new_basepairs.append(newconE)
        
    if (  (con.iunit1 == unit_P and con.iunit2 != unit_A and con.iunit2 != unit_E and con.iunit2 != unit_m) 
        or(con.iunit2 == unit_P and con.iunit1 != unit_A and con.iunit1 != unit_E and con.iunit1 != unit_m)) :
        newconA = BasePair()
        newconE = BasePair()
        if con.iunit1 == unit_P :
            newconA.imp1 = imp_pair_A[0] + (con.imp1 - imp_pair_P[0])
            newconE.imp1 = imp_pair_E[0] + (con.imp1 - imp_pair_P[0])
            newconA.imp2 = con.imp2
            newconE.imp2 = con.imp2
            newconA.iunit1 = unit_A
            newconE.iunit1 = unit_E
            newconA.iunit2 = con.iunit2
            newconE.iunit2 = con.iunit2
        else :
            newconA.imp2 = imp_pair_A[0] + (con.imp2 - imp_pair_P[0])
            newconE.imp2 = imp_pair_E[0] + (con.imp2 - imp_pair_P[0])
            newconA.imp1 = con.imp1
            newconE.imp1 = con.imp1
            newconA.iunit2 = unit_A
            newconE.iunit2 = unit_E
            newconA.iunit1 = con.iunit1
            newconE.iunit1 = con.iunit1
        newconA.imp1un = con.imp1un
        newconE.imp1un = con.imp1un
        newconA.imp2un = con.imp2un
        newconE.imp2un = con.imp2un
        newconA.type = con.type
        newconE.type = con.type
        newconA.coef = con.coef
        newconE.coef = con.coef
        newconA.factor = con.factor
        newconE.factor = con.factor
        newconA.dummy = con.dummy
        newconE.dummy = con.dummy
        newconA.correct_mgo = con.correct_mgo
        newconE.correct_mgo = con.correct_mgo
        newconA.native = con.native
        newconE.native = con.native
        newconA.nhb = con.nhb
        newconE.nhb = con.nhb
        new_basepairs.append(newconA)
        new_basepairs.append(newconE)
        
    if (  (con.iunit1 == unit_E and con.iunit2 != unit_A and con.iunit2 != unit_P and con.iunit2 != unit_m) 
        or(con.iunit2 == unit_E and con.iunit1 != unit_A and con.iunit1 != unit_P and con.iunit1 != unit_m)) :
        newconA = BasePair()
        newconP = BasePair()
        if con.iunit1 == unit_E :
            newconA.imp1 = imp_pair_A[0] + (con.imp1 - imp_pair_E[0])
            newconP.imp1 = imp_pair_P[0] + (con.imp1 - imp_pair_E[0])
            newconA.imp2 = con.imp2
            newconP.imp2 = con.imp2
            newconA.iunit1 = unit_A
            newconP.iunit1 = unit_P
            newconA.iunit2 = con.iunit2
            newconP.iunit2 = con.iunit2
        else :
            newconA.imp2 = imp_pair_A[0] + (con.imp2 - imp_pair_E[0])
            newconP.imp2 = imp_pair_P[0] + (con.imp2 - imp_pair_E[0])
            newconA.imp1 = con.imp1
            newconP.imp1 = con.imp1
            newconA.iunit2 = unit_A
            newconP.iunit2 = unit_P
            newconA.iunit1 = con.iunit1
            newconP.iunit1 = con.iunit1
        newconA.imp1un = con.imp1un
        newconP.imp1un = con.imp1un
        newconA.imp2un = con.imp2un
        newconP.imp2un = con.imp2un
        newconA.type = con.type
        newconP.type = con.type
        newconA.coef = con.coef
        newconP.coef = con.coef
        newconA.factor = con.factor
        newconP.factor = con.factor
        newconA.dummy = con.dummy
        newconP.dummy = con.dummy
        newconA.correct_mgo = con.correct_mgo
        newconP.correct_mgo = con.correct_mgo
        newconA.native = con.native
        newconP.native = con.native
        newconA.nhb = con.nhb
        newconP.nhb = con.nhb
        new_basepairs.append(newconA)
        new_basepairs.append(newconP)

ninfo.contacts.extend(new_contacts)
ninfo.basepairs.extend(new_basepairs)

ninfo.reassign_id()
f_out_ninfo.write_all(ninfo)
f_out_ninfo.close()

imp_set = set()
dupl_set = set()
for con in ninfo.contacts :
    if con.imp1 < con.imp2 :
        imps = (con.imp1, con.imp2)
    else :
        imps = (con.imp2, con.imp1)
    
    if imps in imp_set :
        dupl_set.add(imps)
    else :
        imp_set.add(imps)

for imps in dupl_set:
    print '#################### Warning: contact duplication'
    print '%i %i' % imps
    for con in ninfo.contacts:
        if con.imp1 in imps and con.imp2 in imps :
            con.show()
            print ''

imp_set = set()
dupl_set = set()
for con in ninfo.basepairs :
    if con.imp1 < con.imp2 :
        imps = (con.imp1, con.imp2)
    else :
        imps = (con.imp2, con.imp1)
    
    if imps in imp_set :
        dupl_set.add(imps)
    else :
        imp_set.add(imps)

for imps in dupl_set:
    print '#################### Warning: contact duplication'
    print '%i %i' % imps
    for con in ninfo.basepairs:
        if con.imp1 in imps and con.imp2 in imps :
            con.show()
            print ''

