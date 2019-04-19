#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2012/01/24
@author: Naoto Hori
'''

import sys
from cafysis.file_io.ninfo import NinfoFile
from cafysis.elements.ninfo import NinfoSet

if __name__ == '__main__':
    if len(sys.argv) != 10:
        print('Usage: SCRIPT [ninfo (one-by-one)] [unit] [dist] [side] [unit] [dist] [side] [edc (post_edc.out)] [out]')
        print('')
        sys.exit(2)
        
    fninfo_in = NinfoFile(sys.argv[1])
    fninfo_in.open_to_read()
    ns = NinfoSet()
    fninfo_in.read_all(ns)
    ns.update_info()
    fninfo_in.close()
    
    fninfo_out = NinfoFile(sys.argv[-1])
    fninfo_out.open_to_write()
    
    unit_aln1 = int(sys.argv[2])
    aln1 = {}
    side = int(sys.argv[4])
    for line in open(sys.argv[3]):
        linesp = line.split()
        res1 = int(linesp[0]) # Dunkle res
        mp1 = int(linesp[1]) # Dunkle mp
        res2 = int(linesp[2]) # Jenner res
        mp2 = int(linesp[3]) # Jenner mp
        #linesp[4] # dist
        if side == 1:
            aln1[mp1] = mp2
        elif side == 2:
            aln1[mp2] = mp1
        else:
            print('"side" is not acceptable. See usage')
            sys.exit(2)
    
    unit_aln2 = int(sys.argv[5])
    aln2 = {}
    side = int(sys.argv[7])
    for line in open(sys.argv[6]):
        linesp = line.split()
        res1 = int(linesp[0]) # Dunkle res
        mp1 = int(linesp[1]) # Dunkle mp
        res2 = int(linesp[2]) # Jenner res
        mp2 = int(linesp[3]) # Jenner mp
        #linesp[4] # dist
        #if first:
        #    delta1 = mp1 - 1
        #    delta2 = mp2 - 1
        #    first = False
        if side == 1:
            #aln[mp1 - delta1] = mp2 - delta2
            aln2[mp1] = mp2
        elif side == 2:
            #aln[mp2 - delta2] = mp1 - delta1
            aln2[mp2] = mp1
        else:
            print('"side" is not acceptable. See usage')
            sys.exit(2)
            
    class pairinfo:
        def __init__(self):
            self.ichain1 = 0
            self.ires1 = 0
            self.ires1_u = 0
            self.res1_name = ''
            self.imp1 = 0
            self.imp1_u = 0
            self.mp1_name = ''
            self.ichain2 = 0
            self.ires2 = 0
            self.ires2_u = 0
            self.res2_name = ''
            self.imp2 = 0
            self.imp2_u = 0
            self.mp2_name = ''
            self.iaa1 = 0
            self.iaa2 = 0
            self.eneryg = 0.0
            
    pairlist = {}
    for line in open(sys.argv[8],'r'):
        if len(line) == 1:
            continue
        if line.find('#') != -1:
            continue    
        linesp = line.split()
        pi = pairinfo()
        pi.ichain1 = int(linesp[0])
        pi.ires1 = int(linesp[1])
        pi.ires1_u = int(linesp[2])
        pi.res1_name = linesp[3].strip()
        pi.imp1 = int(linesp[4])
        pi.imp1_u = int(linesp[5])
        pi.mp1_name = linesp[6].strip()
        # linesp[7] is '|'
        pi.ichain2 = int(linesp[8])
        pi.ires2 = int(linesp[9])
        pi.ires2_u = int(linesp[10])
        pi.res2_name = linesp[11].strip()
        pi.imp2 = int(linesp[12])
        pi.imp2_u = int(linesp[13])
        pi.mp2_name = linesp[14].strip()
        # linesp[15] is '|'
        pi.iaa1 = int(linesp[16])
        pi.iaa2 = int(linesp[17])
        # linesp[18] is '|'
        pi.energy = float(linesp[19])
        pairlist[(pi.imp1_u, pi.imp2_u)] = pi
            
    for con in ns.contacts:
        iunit1 = con.iunit1
        iunit2 = con.iunit2
        if unit_aln1 == iunit1 and unit_aln2 == iunit2:
            aln_mp1_u = aln1[con.imp1un]
            aln_mp2_u = aln2[con.imp2un]
        elif unit_aln1 == iunit2 and unit_aln2 == iunit2:
            aln_mp1_u = aln2[con.imp1un]
            aln_mp2_u = aln1[con.imp2un]
        else:
            print('Error:')
            print('unit_aln1:',unit_aln1)
            print('unit_aln2:',unit_aln2)
            print('con.iunit1:',iunit1)
            print('con.iunit2:',iunit2)
            sys.exit(2)
            
        if (aln_mp1_u,aln_mp2_u) in pairlist:
            energy = pairlist[(aln_mp1_u,aln_mp2_u)].energy
        elif (aln_mp2_u,aln_mp1_u) in pairlist:
            energy = pairlist[(aln_mp2_u,aln_mp1_u)].energy
        else:
            energy = 0.0
            print('caution: pair is not found.')
            print('icon=',con.id)
        
        if energy < 0.0:
            con.factor = -1.0 * energy
        else:
            con.factor = 0.0
    
    fninfo_out.write_all(ns)
    fninfo_out.close()
            