#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2012/01/21
@author: Naoto Hori
'''

import sys
from lop.elements.ninfo import NinfoSet
from lop.file_io.ninfo import NinfoFile

if __name__ == '__main__':
    if len(sys.argv) != 6:
        print('Usage: %SCRIPT [ninfo] [target unit] [aln file] [side] [OUT ninfo]')
        print('       side=1  old -> new')
        print('       side=2  new -> old')
        sys.exit(2)

    target = int(sys.argv[2])
    side = int(sys.argv[4])
    
    aln = {}
    #first = True
    for line in open(sys.argv[3], 'r'):
        linesp = line.split()
        res1 = int(linesp[0]) # Dunkle res # Hantai???
        mp1 = int(linesp[1]) # Dunkle mp
        res2 = int(linesp[2]) # Jenner res # Hantai???
        mp2 = int(linesp[3]) # Jenner mp
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
            print('"side" is not acceptable. See usage')
            sys.exit(2)
            
    ns = NinfoSet()
    file_ninfo = NinfoFile(sys.argv[1])
    file_ninfo.open_to_read()
    file_ninfo.read_all(ns)
    file_ninfo.close()
    
    ns_new = NinfoSet()
    
    # BondLengths
    add = False
    for ni in ns.bondlengths:
        if ni.iunit1 == target:
            if ni.imp1un in aln:
                ni.imp1un = aln[ni.imp1un]
                add = True
        if ni.iunit2 == target:
            if ni.imp2un in aln:
                ni.imp2un = aln[ni.imp2un]
                add = True
        if add:
            ns_new.bondlengths.append(ni)
            add = False
            
    # BondAngles
    add = False
    for ni in ns.bondangles:
        if ni.iunit1 == target:
            if ni.imp1un in aln:
                ni.imp1un = aln[ni.imp1un]
                add = True
        if ni.iunit2 == target:
            if ni.imp2un in aln:
                ni.imp2un = aln[ni.imp2un]
                add = True
        if add:
            ns_new.bondangles.append(ni)
            add = False

    # Dihedrals
    add = False
    for ni in ns.dihedrals:
        if ni.iunit1 == target:
            if ni.imp1un in aln:
                ni.imp1un = aln[ni.imp1un]
                add = True
        if ni.iunit2 == target:
            if ni.imp2un in aln:
                ni.imp2un = aln[ni.imp2un]
                add = True
        if add:
            ns_new.dihedrals.append(ni)
            add = False

    # Contacts
    add = False
    for ni in ns.contacts:
        if ni.iunit1 == target:
            if ni.imp1un in aln:
                ni.imp1un = aln[ni.imp1un]
                add = True
        if ni.iunit2 == target:
            if ni.imp2un in aln:
                ni.imp2un = aln[ni.imp2un]
                add = True
        if add:
            ns_new.contacts.append(ni)
            add = False

    # BasePairs
    add = False
    for ni in ns.basepairs:
        if ni.iunit1 == target:
            if ni.imp1un in aln:
                ni.imp1un = aln[ni.imp1un]
                add = True
        if ni.iunit2 == target:
            if ni.imp2un in aln:
                ni.imp2un = aln[ni.imp2un]
                add = True
        if add:
            ns_new.basepairs.append(ni)
            add = False

    # BaseStacks
    add = False
    for ni in ns.basestacks:
        if ni.iunit1 == target:
            if ni.imp1un in aln:
                ni.imp1un = aln[ni.imp1un]
                add = True
        if ni.iunit2 == target:
            if ni.imp2un in aln:
                ni.imp2un = aln[ni.imp2un]
                add = True
        if add:
            ns_new.basestacks.append(ni)
            add = False
        
    out_ninfo = NinfoFile(sys.argv[5])
    out_ninfo.open_to_write()
    out_ninfo.write_all(ns_new)
    out_ninfo.close()