#!/usr/bin/env python
#vim:fileencoding=UTF-8
'''
Modified based on dcd_body_origin_PBC.py on 2017/12/16

Created on 2017/10/13
@author: Naoto Hori

Translate all coordinates so that the entire molecule is inside the box.
To do so, firstly find most lateral coordinates.

Particles outside the box will be wrapped into the box in the origin.
'''

from cafysis.mtx_coord_transform import mtx_crd_transform
from cafysis.file_io.dcd import DcdFile
from cafysis.file_io.pdb import PdbFile
import py_bestfit
import numpy as np
from numpy import zeros, float64
import sys
import math
import copy

## Threshold of distance of neighboring beads
MAXD = 50.0

if len(sys.argv) != 7:
    print 'Usage: SCRIPT [input DCD] [input PDB] [ID domain begin] [ID domain end] [Box size] [output DCD]'
    sys.exit(2)

ID_DOM_INI = int(sys.argv[3]) - 1  # 重心を求める際に必要
ID_DOM_END = int(sys.argv[4]) - 1
BOXSIZE = float(sys.argv[5])

BOXMAX = 0.5 * BOXSIZE
BOXMIN = -0.5 * BOXSIZE


### Files
dcd = DcdFile(sys.argv[1])
dcd.open_to_read()

pdb = PdbFile(sys.argv[2])
pdb.open_to_read()

dcd_out = DcdFile(sys.argv[-1])
dcd_out.open_to_write()


### Read the reference pdb
ref_chains = pdb.read_all()

num_atom = 0
for chain in ref_chains :
    num_atom += chain.num_atom()

ref = zeros((3, num_atom), dtype=float64, order='F')
i = 0
for chain in ref_chains :
    for residue in chain.residues:
        for atom in residue.atoms :
            (ref[0][i], ref[1][i], ref[2][i]) = atom.xyz.get_as_tuple()
            i += 1
pdb.close()

ref_idx = [i for i in range(ID_DOM_INI, ID_DOM_END+1)]
pre_idx = [i for i in range(ID_DOM_INI, ID_DOM_END+1)]


# header
dcd.read_header()
header = dcd.get_header()
dcd_out.set_header(header)
dcd_out.write_header()

nmp = header.nmp_real

def find_max_min_PBC(d):

    max_xyz = np.copy(d[0,:])
    min_xyz = np.copy(d[0,:])
    pre     = np.copy(d[0][:])
    add     = np.zeros( (3,))

    nmp, _ = d.shape

    for ixyz in range(nmp):

        xyz = d[ixyz,0:3]

        for i in range(3):
            x = xyz[i] + add[i]
            if x - pre[i] > MAXD:
                x      += - BOXSIZE
                add[i] += - BOXSIZE
            elif x - pre[i] < -MAXD:
                x      +=   BOXSIZE
                add[i] +=   BOXSIZE
            pre[i] = x

            if max_xyz[i] < x:
                max_xyz[i] = x
            if x < min_xyz[i]:
                min_xyz[i] = x

    L = max_xyz - min_xyz
    # L.shape is (3,)

    return max_xyz, min_xyz, L

def wrap(d):
    for i in range(len(d)):
        for j in range(3):
            p = d[i][j]
            if p > BOXMAX:
                d[i][j] = p - BOXSIZE * (int((p - BOXMAX)/BOXSIZE) + 1)
            elif p < BOXMIN:
                d[i][j] = p + BOXSIZE * (int((BOXMIN - p)/BOXSIZE) + 1)

iframe = 0
while dcd.has_more_data() :

    data = dcd.read_onestep_np()
    
    ##########################################################
    ### Wrapping the molecule to do bestfit

    # Find most lateral coordinates
    max_xyz, min_xyz, L = find_max_min_PBC(data[ID_DOM_INI:ID_DOM_END+1])

    if L[0] > BOXSIZE or L[1] > BOXSIZE or L[2] > BOXSIZE:
        print ('Warning: L exceeds BOXSIZE at frame %i' % iframe)

    mtx_move = mtx_crd_transform()
    #mtx_move.reset()
    mtx_move.translation(0.5*L[0] - max_xyz[0],
                         0.5*L[1] - max_xyz[1],
                         0.5*L[2] - max_xyz[2])
    
    mtx_move.do_to_ndarray(data)


    ##########################################################
    ### Align to the reference PDB
    (post, rmsd, ier, rot, center_ref, center_pre) = py_bestfit.bestfit(ref, data.T,
                                                                        ref_idx, pre_idx)
    data = post.T


    ##########################################################
    ### Wrapping again because some part may be outside because of bestfit

    # Find most lateral coordinates
    max_xyz, min_xyz, L = find_max_min_PBC(data[ID_DOM_INI:ID_DOM_END+1])

    if L[0] > BOXSIZE or L[1] > BOXSIZE or L[2] > BOXSIZE:
        print ('Warning: L exceeds BOXSIZE at frame %i' % iframe)

    mtx_move = mtx_crd_transform()
    #mtx_move.reset()
    mtx_move.translation(0.5*L[0] - max_xyz[0],
                         0.5*L[1] - max_xyz[1],
                         0.5*L[2] - max_xyz[2])
    
    mtx_move.do_to_ndarray(data)

    ##########################################################
    # Wrap particles outside the box at the orgin
    wrap(data)

    dcd_out.write_onestep(data) 

    iframe += 1
    
dcd.close()
dcd_out.close()
