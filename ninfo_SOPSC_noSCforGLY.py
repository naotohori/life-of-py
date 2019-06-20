#!/usr/bin/env python

'''
Created on 2017/11/15
@author: Naoto Hori
'''

FENE_R0 = 2.0
FENE_K = 20.0

Rc = 8.0

EPS_B_B = 0.45
EPS_B_S = 0.45
COEF_S_S = 0.3  # eps = COEF_S_S * abs(BT - 0.7)

SEP_B_B = 3  # residue (i) and residue (i+3)
SEP_B_S = 2  # residue (i) and residue (i+2)
SEP_S_S = 2  # residue (i) and residue (i+2)

import sys
from cafysis.file_io.pdb import PdbFile
from cafysis.file_io.ninfo import NinfoFile
from cafysis.elements.ninfo import NinfoSet, Fene, LJ
from cafysis.para.BetancourtThirumalai import BTmatrix

if len(sys.argv) != 3:
    print('Usage: SCRIPT [input PDB] [output ninfo]')
    sys.exit(2)

p = PdbFile(sys.argv[1])
p.open_to_read()

chains = p.read_all()

ns = NinfoSet()


a_BB_pre = None

for ic, c in enumerate(chains):
    for r in c.residues:

        num_a = len(r.atoms)
        a_BB = None
        a_SC = None

        for a in r.atoms:
            if a.res_name.strip() == 'GLY':
                if num_a != 1:
                    print('Error: num_a != 1')
                    print(a.name, a.res_name, a.res_seq)
                    sys.exit(2)
            else:
                if num_a != 2:
                    print('Error: num_a != 2')
                    print(a.name, a.res_name, a.res_seq)
                    sys.exit(2)

            if a.name.strip() == 'B':
                a_BB = a
            elif a.name.strip() == 'S':
                a_SC = a
            else:
                print('Error: unknown particle type')
                print(a.name, a.res_name, a.res_seq)
                sys.exit(2)

        ## bond between B(i-1) and B(i)
        if a_BB_pre is not None:
            f = Fene(id=None,iunit1=ic+1,iunit2=ic+1, imp1=a_BB_pre.serial, imp2=a_BB.serial,
                     imp1un=a_BB_pre.serial,imp2un=a_BB.serial,
                     native=a_BB.xyz.distance(a_BB_pre.xyz),
                     factor = FENE_R0, coef = FENE_K,
                     correct_mgo=None, type_str=None)
            ns.fenes.append(f)
        a_BB_pre = a_BB

        ## bond between B(i) and S(i)
        if num_a == 2:
            if a_SC is None:
                print('Error: a_SC is None')
                sys.exit(2)
            ni = Fene(id=None,iunit1=ic+1,iunit2=ic+1, imp1=a_BB.serial, imp2=a_SC.serial,
                     imp1un=a_BB.serial,imp2un=a_SC.serial,
                     native=a_BB.xyz.distance(a_SC.xyz),
                     factor = FENE_R0, coef = FENE_K,
                     correct_mgo=None, type_str=None)
            ns.fenes.append(ni)

'''
Suppose only 1 chain
'''
num_atom = chains[0].num_atom()

for ia1 in range(num_atom):
    
    a1 = chains[0].get_atom(ia1)
    res1 = a1.res_seq
    name1 = a1.name.strip()

    for ia2 in range(ia1+1, num_atom):

        a2 = chains[0].get_atom(ia2)
        res2 = a2.res_seq
        name2 = a2.name.strip()

        # Backbone - Backbone 
        if name1 == 'B' and name2 == 'B':
            if abs(res2-res1) >= SEP_B_B:
                dist = a1.xyz.distance(a2.xyz)
                if dist <= Rc:
                    ni = LJ(id=None,iunit1=1,iunit2=1,
                            imp1=a1.serial, imp2=a2.serial,
                            imp1un=a1.serial, imp2un=a2.serial,
                            native = dist, coef = EPS_B_B,
                            factor=None, correct_mgo=None, type_str=None)
                    ns.LJs.append(ni)

        # Sidechain - Sidechain 
        elif name1 == 'S' and name2 == 'S':
            if abs(res2-res1) >= SEP_S_S:
                dist = a1.xyz.distance(a2.xyz)
                if dist <= Rc:

                    if (a1.res_name, a2.res_name) in BTmatrix:
                        BTeps = BTmatrix[(a1.res_name, a2.res_name)]
                    elif (a2.res_name, a1.res_name) in BTmatrix:
                        BTeps = BTmatrix[(a2.res_name, a1.res_name)]
                    else:
                        print('Erorr: either a1.res_name=%s or a2.res_name=%s is wrong' % (a1.res_name, a2.res_name))
                        sys.exit(2)

                    eps = COEF_S_S * abs(BTeps - 0.7)

                    ni = LJ(id=None,iunit1=1,iunit2=1,
                            imp1=a1.serial, imp2=a2.serial,
                            imp1un=a1.serial, imp2un=a2.serial,
                            native = dist, coef = eps,
                            factor=None, correct_mgo=None, type_str=None)
                    ns.LJs.append(ni)

        # Backbone - Sidechain
        elif (name1 == 'B' and name2 == 'S') or (name1 == 'S' and name2 == 'B'):
            if abs(res2-res1) >= SEP_B_S:
                dist = a1.xyz.distance(a2.xyz)
                if dist <= Rc:
                    ni = LJ(id=None,iunit1=1,iunit2=1,
                            imp1=a1.serial, imp2=a2.serial,
                            imp1un=a1.serial, imp2un=a2.serial,
                            native = dist, coef = EPS_B_S,
                            factor=None, correct_mgo=None, type_str=None)
                    ns.LJs.append(ni)
        else:
            print("Error: either name1=%s or name2=%s is unknown" % (name1,name2))
            sys.exit(2)

nf = NinfoFile(sys.argv[-1])
nf.open_to_write()
nf.write_all(ns)
nf.close()
