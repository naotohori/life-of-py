#!/usr/bin/env python

import sys
from file_pdb import PdbFile

if len(sys.argv) != 5 :
    print ('\n Usage: SCRIPT [list file] [input PDB1] [input PDB2] [output]\n')
    sys.exit(2)
    
f_list = open(sys.argv[1],'r')

f_pdb = PdbFile(sys.argv[2])
f_pdb.open_to_read()
chains1 = f_pdb.read_all()
f_pdb.close

f_pdb = PdbFile(sys.argv[3])
f_pdb.open_to_read()
chains2 = f_pdb.read_all()
f_pdb.close

f_out = open(sys.argv[4], 'w')

if len(chains1) != len(chains2) :
    print(("Error: len(chains1)(=%i)  !=  len(chains2)(=%i)" %(len(chains1), len(chains2))))
    sys.exit(2)
# !!! current version is for only single chain !!!
if len(chains1) != 1 or len(chains2) != 1:
    print ("Error: len(chains1) != 1 or len(chains2) != 1" )
    sys.exit(2)
    
res1_to_res2 = {}
for line in f_list :
    if line.find('#') != -1 :
        continue
    linesp = line.split()
    ires1 = int(linesp[0])
    ires2 = int(linesp[1])
    res1_to_res2[ires1] = ires2
    

c1 = chains1[0]
c2 = chains2[0]

data = ()
for ir,r1 in enumerate(c1.residues):
    ires1 = ir + 1
    if not ires1 in res1_to_res2 :
        continue
    ires2 = res1_to_res2[ires1]
    r2 = c2.residues[ires2-1]
    
    for a1 in r1.atoms:
        for a2 in r2.atoms :
            if a1.name in (' S  ', ' P  ') or a2.name in (' S  ', ' P  ') :
                if a1.name == a2.name :
                    distance = a1.xyz.distance(a2.xyz)
                    f_out.write('%5i %6i %5i %6i %6.2f\n' % (ires1, a1.serial, ires2, a2.serial, distance))
            else :
                distance = a1.xyz.distance(a2.xyz)
                f_out.write('%5i %6i %5i %6i %6.2f\n' % (ires1, a1.serial, ires2, a2.serial, distance))

f_out.close()