#!/usr/bin/env python
#vim:fileencoding=UTF-8
'''
Created on 2024/12/12
@author: Huong
adopted from pdb_centeroid_origin_PBC.py from Naoto

Unwrap particles outside the box.

Note: Box size need to be bigger than individual chain sizes.

'''

from lop.mtx_coord_transform import mtx_crd_transform
from lop.file_io.pdb import PdbFile
import sys
import math
import copy

if len(sys.argv) != 4:
    print('Usage: SCRIPT [input PDB] [Box size] [output PDB]')
    sys.exit(2)

BOXSIZE = float(sys.argv[2])

pdb = PdbFile(sys.argv[1])
pdb.open_to_read()
pdb_out = PdbFile(sys.argv[-1])
pdb_out.open_to_write()

# header
chains = pdb.read_all()
data = []
for c in chains:
    for r in c.residues:
        for a in r.atoms:
            data.append([a.xyz.x,a.xyz.y,a.xyz.z])


def calc_com_PBC(d):
    '''
    See https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions
    '''
    cos = [0.0] * 3
    sin = [0.0] * 3
    for xyz in d:
        for i in range(3):
            theta = (xyz[i] ) / BOXSIZE * 2 * math.pi
            cos[i] += math.cos(theta)
            sin[i] += math.sin(theta)

    com = [0.0] * 3
    for i in range(3):
        cos[i] = cos[i] / float(len(d))
        sin[i] = sin[i] / float(len(d))
        theta = math.atan2(-sin[i],-cos[i]) + math.pi
        com[i] = 0.5 * BOXSIZE * theta / math.pi 

    return com

def wrap(d):
    for i in range(len(d)):
        for j in range(3):
            p = d[i][j]
            if p > BOXMAX[j]:
                d[i][j] = p - BOXSIZE * (int((p - BOXMAX[j])/BOXSIZE) + 1)
            elif p < BOXMIN[j]:
                d[i][j] = p + BOXSIZE * (int((BOXMIN[j] - p)/BOXSIZE) + 1)

###########################################################
nmp = 0
i = 0
for c in chains:
    ID_DOM_INI = nmp
    nmp += c.num_atom()
    ID_DOM_END = nmp - 1
#    print(str(ID_DOM_INI)+'_'+str(ID_DOM_END))
    COM = calc_com_PBC(data[ID_DOM_INI:ID_DOM_END+1])

    ####################################3
    # Set box boundaries
    BOXMAX = [0.0] * 3
    BOXMIN = [0.0] * 3
    for j in range(3):
        BOXMAX[j] = COM[j] + 0.5 * BOXSIZE
        BOXMIN[j] = COM[j] - 0.5 * BOXSIZE

    ##########################################################
    # Wrap particles outside the box centered at COM
    wrap(data)

    for r in c.residues:
        for a in r.atoms:
            a.xyz.x = data[i][0]
            a.xyz.y = data[i][1]
            a.xyz.z = data[i][2]
            i += 1
    
pdb_out.write_all(chains)
pdb_out.close()
