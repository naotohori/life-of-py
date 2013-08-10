#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2012/01/21
@author: Naoto Hori
'''

import sys
from cafysis.elements.ninfo import NinfoSet
from cafysis.file_io.ninfo import NinfoFile

if __name__ == '__main__':
    if len(sys.argv) != 6:
        print 'Usage: %SCRIPT [ninfo] [target unit] [aln file] [side] [OUT ninfo]'
        print '       side=1  old -> new'
        print '       side=2  new -> old'
        sys.exit(2)

    target = int(sys.argv[2])
    side = int(sys.argv[4])
    
    aln = {}
    #first = True
    for line in open(sys.argv[3], 'r'):
        linesp = line.split()
        #res1 = int(linesp[0])
        mp1 = int(linesp[0])
        #res2 = int(linesp[2])
        mp2 = int(linesp[1])
        #linesp[4] # dist
        #if first:
        #    delta1 = mp1 - 1
        #    delta2 = mp2 - 1
        #    first = False
        if side == 1:
            #aln[mp1 - delta1] = mp2 - delta2
            aln[mp1] = mp2
        elif side == 2:
            #aln[mp2 - delta2] = mp1 - delta1
            aln[mp2] = mp1
        else:
            print '"side" is not acceptable. See usage'
            sys.exit(2)
            
    ns = NinfoSet()
    file_ninfo = NinfoFile(sys.argv[1])
    file_ninfo.open_to_read()
    file_ninfo.read_all(ns)
    file_ninfo.close()
    
    ns_new = NinfoSet()
    
    
    # 両方のユニットともターゲットのときは、
    # 両方のIDがalnに含まれていないとaddしない
    
    # BondLengths
    add1 = False
    add2 = False
    for ni in ns.bondlengths:
        if ni.iunit1 == target:
            if ni.imp1un in aln:
                ni.imp1un = aln[ni.imp1un]
                add1 = True
            if ni.imp2un in aln:
                ni.imp2un = aln[ni.imp2un]
                add2 = True
        if ni.iunit1 == target:
            if add1 and add2:
                ns_new.bondlengths.append(ni)
        else:
            ns_new.bondlengths.append(ni)
        add1 = False
        add2 = False
            
    # BondAngles
    add1 = False
    add2 = False
    add3 = False
    for ni in ns.bondangles:
        if ni.iunit1 == target:
            if ni.imp1un in aln:
                ni.imp1un = aln[ni.imp1un]
                add1 = True
            if ni.imp2un in aln:
                ni.imp2un = aln[ni.imp2un]
                add2 = True
            if ni.imp3un in aln:
                ni.imp3un = aln[ni.imp3un]
                add3 = True
        if ni.iunit1 == target:
            if add1 and add2 and add3:
                ns_new.bondangles.append(ni)
        else:
            ns_new.bondangles.append(ni)
        add1 = False
        add2 = False
        add3 = False
        
    # AICG13
    add1 = False
    add2 = False
    add3 = False
    for ni in ns.aicg13s:
        if ni.iunit1 == target:
            if ni.imp1un in aln:
                ni.imp1un = aln[ni.imp1un]
                add1 = True
            if ni.imp2un in aln:
                ni.imp2un = aln[ni.imp2un]
                add2 = True
            if ni.imp3un in aln:
                ni.imp3un = aln[ni.imp3un]
                add3 = True
        if ni.iunit1 == target:
            if add1 and add2 and add3:
                ns_new.aicg13s.append(ni)
        else:
            ns_new.aicg13s.append(ni)
        add1 = False
        add2 = False
        add3 = False

    # Dihedrals
    add1 = False
    add2 = False
    add3 = False
    add4 = False
    for ni in ns.dihedrals:
        if ni.iunit1 == target:
            if ni.imp1un in aln:
                ni.imp1un = aln[ni.imp1un]
                add1 = True
            if ni.imp2un in aln:
                ni.imp2un = aln[ni.imp2un]
                add2 = True
            if ni.imp3un in aln:
                ni.imp3un = aln[ni.imp3un]
                add3 = True
            if ni.imp4un in aln:
                ni.imp4un = aln[ni.imp4un]
                add4 = True
        if ni.iunit1 == target:
            if add1 and add2 and add3 and add4:
                ns_new.dihedrals.append(ni)
        else:
            ns_new.dihedrals.append(ni)
        add1 = False
        add2 = False
        add3 = False
        add4 = False

    # AicgDih
    add1 = False
    add2 = False
    add3 = False
    add4 = False
    for ni in ns.aicgdihs:
        if ni.iunit1 == target:
            if ni.imp1un in aln:
                ni.imp1un = aln[ni.imp1un]
                add1 = True
            if ni.imp2un in aln:
                ni.imp2un = aln[ni.imp2un]
                add2 = True
            if ni.imp3un in aln:
                ni.imp3un = aln[ni.imp3un]
                add3 = True
            if ni.imp4un in aln:
                ni.imp4un = aln[ni.imp4un]
                add4 = True
        if ni.iunit1 == target:
            if add1 and add2 and add3 and add4:
                ns_new.aicgdihs.append(ni)
        else:
            ns_new.aicgdihs.append(ni)
        add1 = False
        add2 = False
        add3 = False
        add4 = False

    # Contacts
    add1 = False
    add2 = False
    for ni in ns.contacts:
        if ni.iunit1 == target:
            if ni.imp1un in aln:
                ni.imp1un = aln[ni.imp1un]
                add1 = True
        if ni.iunit2 == target:
            if ni.imp2un in aln:
                ni.imp2un = aln[ni.imp2un]
                add2 = True
        if ((ni.iunit1 == target and not add1) or
            (ni.iunit2 == target and not add2)):
            pass
        else:
            ns_new.contacts.append(ni)
        add1 = False
        add2 = False

    # BasePairs
    add1 = False
    add2 = False
    for ni in ns.basepairs:
        if ni.iunit1 == target:
            if ni.imp1un in aln:
                ni.imp1un = aln[ni.imp1un]
                add1 = True
        if ni.iunit2 == target:
            if ni.imp2un in aln:
                ni.imp2un = aln[ni.imp2un]
                add2 = True
        if ((ni.iunit1 == target and not add1) or
            (ni.iunit2 == target and not add2)):
            pass
        else:
            ns_new.basepairs.append(ni)
        add1 = False
        add2 = False

    # BaseStacks
    add1 = False
    add2 = False
    for ni in ns.basestacks:
        if ni.iunit1 == target:
            if ni.imp1un in aln:
                ni.imp1un = aln[ni.imp1un]
                add1 = True
        if ni.iunit2 == target:
            if ni.imp2un in aln:
                ni.imp2un = aln[ni.imp2un]
                add2 = True
        if ((ni.iunit1 == target and not add1) or
            (ni.iunit2 == target and not add2)):
            pass
        else:
            ns_new.basestacks.append(ni)
        add1 = False
        add2 = False
        
    out_ninfo = NinfoFile(sys.argv[5])
    out_ninfo.open_to_write()
    out_ninfo.write_all(ns_new)
    out_ninfo.close()