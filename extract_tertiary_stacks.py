#!/usr/bin/env python

import sys

if len(sys.argv) != 4:
    print("Usage: SCRIPT [input hb.native.dat] [input pdb_detect_basestack.out] [output stack.dat]")
    sys.exit(2)

filename_hb = sys.argv[1]
filename_stack = sys.argv[2]
filename_out = sys.argv[3]

canonical_base_pairs = []
for l in open(filename_hb):
    if l[0:6] == 'NATIVE':
        continue
    if l[0:3] == 'NAT':
        lsp = l.split()
        if len(lsp) < 5:
            print("Error: len(lsp) < 5, this is not expected")
            print("line: %s" % l)
            sys.exit(2)

        i = int(lsp[1])
        j = int(lsp[2])

        if i < j:
            canonical_base_pairs.append((i,j))
        else:
            canonical_base_pairs.append((j,i))

#print(canonical_base_pairs)

f_out = open(filename_out,"w")

for l in open(filename_stack):
    if l[0] == 'S':
        continue
    if l[0] != 'T':
        print("Error: unknown format in stack file")
        print("line: %s" % l)
        sys.exit(2)

    lsp = l.split()
    if len(lsp) != 7:
        print("Error: len(lsp) != 7 in stack file, this is not expected")
        print("line: %s" % l)
        sys.exit(2)

    signed_i = int(lsp[1])
    signed_j = int(lsp[2])

    i = abs(signed_i)
    j = abs(signed_j)
    #print(i,j)

    if j < i:
        i_tmp = i
        i = j
        j = i_tmp

    if ((i-1,j) in canonical_base_pairs) and ((i,j-1) in canonical_base_pairs):
        continue
    elif ((i,j+1) in canonical_base_pairs) and ((i+1,j) in canonical_base_pairs):
        continue
    else:
        f_out.write("{:+d} {:+d}\n".format(signed_i, signed_j))
        

f_out.close()
