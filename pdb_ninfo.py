#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2012/05/01
@author: Naoto Hori
'''

if __name__ == '__main__':
    pass

    
import sys
from cafysis.elements.ninfo import NinfoSet
from cafysis.file_io.ninfo import NinfoFile
from cafysis.file_io.pdb import PdbFile

if len(sys.argv) != 4:
    print ("Usage: SCRIPT [ninfo file] [pdb file] [output file]")
    sys.exit(2)

file_ninfo = NinfoFile(sys.argv[1])
file_ninfo.open_to_read()
ns = NinfoSet()
file_ninfo.read_all(ns)
ns.update_info()
file_ninfo.close()

file_pdb = PdbFile(sys.argv[2])
file_pdb.open_to_read()
chains = file_pdb.read_all()

xyz = []
for c in chains:
    for i in range(c.num_atom()):
        xyz.append(c.get_atom(i).xyz)
    
print(("## Confirmation: number of atoms = %i" % len(xyz)))
        
f_out = open(sys.argv[-1],'w')

f_out.write('#Contacts\n')
for con in ns.contacts:
    conid = con.id
    imp1 = con.imp1
    imp2 = con.imp2
    nat = con.native
    dist = xyz[con.imp1-1].distance(xyz[con.imp2-1])
    f_out.write('%10i %5i %5i %8.3f %8.3f\n' % (conid, imp1,imp2,nat,dist))
    
f_out.write('\n\n#Basepairs\n')
for bp in ns.basepairs:
    bpid = bp.id
    imp1 = bp.imp1
    imp2 = bp.imp2
    nat = bp.native
    dist = xyz[bp.imp1-1].distance(xyz[bp.imp2-1])
    f_out.write('%10i %5i %5i %8.3f %8.3f\n' % (bpid, imp1,imp2,nat,dist))
    
f_out.write('\n\n#Basestacks\n')
for bs in ns.basestacks:
    bsid = bs.id
    imp1 = bs.imp1
    imp2 = bs.imp2
    nat = bs.native
    dist = xyz[bs.imp1-1].distance(xyz[bs.imp2-1])
    f_out.write('%10i %5i %5i %8.3f %8.3f\n' % (bsid, imp1,imp2,nat,dist))
    
f_out.close()
    
    