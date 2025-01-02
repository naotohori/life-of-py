#!/usr/bin/env python

import sys
import string

if len(sys.argv) != 3:
    print("Usage: SCRIPT [input PDB] [output PDB]")
    sys.exit(2)

input_PDB = sys.argv[1]
output_PDB = sys.argv[2]

chain_ids = {}
# chains 1 to 10 to be "0" to "9"
for i in range(10):
    chain_ids[i+1] = f"{i}"

for i in range(11,11+26):
    chain_ids[i] = f"{string.ascii_lowercase[i-11]}"

for i in range(11+26,11+26+26):
    chain_ids[i] = f"{string.ascii_uppercase[i-37]}"

chain_ids[63] = '_'
chain_ids[64] = '-'

#for i in range(1,65):
#    print(i, chain_ids[i])
#sys.exit(0)

f_out = open(output_PDB,'w')

new_chain = True
chainid = 1
for l in open(input_PDB):

    if l.startswith('HETATM') or l.startswith('ATOM  '):
        f_out.write(l[0:20] + ' ' + chain_ids[chainid] + l[22:])

    elif l.startswith('TER'):
        try:
            chainid = int(l[20:22])
        except:
            f_out.write(l)
        else:
            f_out.write(l[0:20] + ' ' + chain_ids[chainid] + l[22:])
        chainid += 1

    else:
        f_out.write(l)
    

f_out.close()
