#!/usr/bin/env python

import sys

"""
Take coordinates from the input xyz file,
and insert the coordinates into the given cif file as a template,
then write them to a new cif file.
"""

if len(sys.argv) != 4:
    print('Usage: SCRIPT (xyz file) (cif file) (out file)')
    sys.exit(2)

nbeads = None
xyzs = []
elements = []
for il, l in enumerate(open(sys.argv[1])):
    if il == 0:
        nbeads = int(l)
    elif il == 1:
        continue
    else:
        lsp = l.split()
        elements.append(lsp[0])
        x = float(lsp[1])
        y = float(lsp[2])
        z = float(lsp[3])
        xyzs.append((x,y,z))

fout = open(sys.argv[3], 'w')

ibead = 0
for l in open(sys.argv[2]):
    lsp = l.strip().split()
    if lsp[0].isdigit():
        for i in range(len(lsp)):
            if i == 8:
                fout.write(f'{xyzs[ibead][0]:.3f} ')
            elif i == 9:
                fout.write(f'{xyzs[ibead][1]:.3f} ')
            elif i == 10:
                fout.write(f'{xyzs[ibead][2]:.3f} ')
            else:
                fout.write(lsp[i] + ' ')
        fout.write('\n')
        assert lsp[-1] == elements[ibead]
        ibead += 1
    else:
        fout.write(l)
fout.close()

assert ibead == nbeads

