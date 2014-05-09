#!/usr/bin/env python

from file_pdb import PdbFile
import sys

if len(sys.argv) != 2 :
    print 'Usage: SCRIPT [input PDB]'
    sys.exit(2)

pdb_cafemol = PdbFile(sys.argv[1])
pdb_cafemol.open_to_read()
cafe_chains = pdb_cafemol.read_all()
pdb_cafemol.close()

#last_res_nums = []
#sum = 0
#for c in cafe_chains:
#    last_res_num = c.get_atom(-1).res_seq
#    sum += last_res_num
#    last_res_nums.append(last_res_num)
#print last_res_nums
#print sum

offset = 0
print('#ic   ires  res_c res    imp name')
for (ic, c) in enumerate(cafe_chains):
    for (ir, r) in enumerate(c.residues):
        for (ia, a) in enumerate(r.atoms):
            print('%3i %6i %6i %s %6i %s'
               % (ic+1, a.res_seq+offset, a.res_seq, a.res_name, a.res_seq, a.name))
    offset += c.get_atom(-1).res_seq
