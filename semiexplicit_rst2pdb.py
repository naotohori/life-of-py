#!/usr/bin/env python

import sys

if len(sys.argv) != 4:
    print("Usage: SCRIPT [rst file] [pdb file] [output]")
    sys.exit(2)
rstfile = sys.argv[1]
pdbfile = sys.argv[2]
outfile = sys.argv[3]

coord = []

for l in open(rstfile):
    if l[0:3] == 'ATM':
        nmp = int(l.split()[1])
        print('ATM: %i' % nmp)
        continue

    if l[0:3] == 'BOX':
        break

    lsp = l.split()
    if len(lsp) != 6:
        print("Error: Unkown format in rst file")
        sys.exit(2)

    coord.append( [float(lsp[0]), float(lsp[1]), float(lsp[2])] )

print('Number of lines in rst data: %i' % len(coord))

if nmp != len(coord):
    print("Error: inconsistent number of lines")
    sys.exit(2)


f_out = open(outfile,"w")

for il, l in enumerate(open(pdbfile)):
    if l[0:4] == 'ATOM':
        f_out.write(l[0:30])
        f_out.write("%8.3f%8.3f%8.3f\n" % tuple(coord[il]))
    else:
        f_out.write(l)

    
f_out.close()

