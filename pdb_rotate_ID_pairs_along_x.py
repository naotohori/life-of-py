#!/usr/bin/env python
# authors: Huong T VU, Naoto Hori 
# Jan 2 2024
# rotating given PDB so given pairs of beads align along vetor x

import sys
import re
import os 
import shutil
import operator
import subprocess
from lop.file_io.pdb import PdbFile
from lop.elements.pdb import Chain, Residue, Atom
import numpy as np
def ranges(i):
    for a, b in itertools.groupby(enumerate(i), lambda pair: pair[1] - pair[0]):
        b = list(b)
        yield b[0][1], b[-1][1]

def get_COM(vec_list):
    c = [sum(i)/len(vec_list) for i in zip(*vec_list)]
    return c

def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix

# vec1 = [2, 3, 2.5]
# vec2 = [-3, 1, -3.4]

# mat = rotation_matrix_from_vectors(vec1, vec2)
# vec1_rot = mat.dot(vec1)
# assert np.allclose(vec1_rot/np.linalg.norm(vec1_rot), vec2/np.linalg.norm(vec2))

def rotation_matrix_align_to_x_axis(vec):
    """ Find the rotation matrix that aligns vec to x axis
    :param vec: A 3d "source" vector
    :return mat: A transform matrix (3x3) which when applied to vec, aligns it with [1.0, 0.0, 0.0].
    """ 
    vecx = [1.0, 0.0, 0.0]
    rotation_matrix = rotation_matrix_from_vectors(vec, vecx)
    return rotation_matrix

pdbpath = sys.argv[1]
ID_pairs_first = [int(i) for i in sys.argv[2:-1:2]]
ID_pairs_last = [int(i) for i in sys.argv[3:-1:2]]
if len(ID_pairs_first)==0:
    print ('\n Usage: SCRIPT [PDB in] [ID pair 1 first] [ID pair 1 last] ... [ID pair i last] [PDB out]')
    sys.exit(2)
if len(ID_pairs_first) != len(ID_pairs_last):
    print("ID_pairs list has to have even argvs (pairs)")
    sys.exit(2)
outfilename = sys.argv[-1]
    
f_pdb = PdbFile(pdbpath)
print('Aligning '+pdbpath+' ...')
f_pdb.open_to_read()
chains = f_pdb.read_all_and_close()

xyz_pairs_first = []
xyz_pairs_last = []
for c in chains :
    for r in c.residues :
        for a in r.atoms :
            if a.serial in ID_pairs_first:
                xyz_pairs_first.append(a.xyz.get_as_list())
            elif a.serial in ID_pairs_last:
                xyz_pairs_last.append(a.xyz.get_as_list())
vec1 = [ai - bi for ai, bi in zip(get_COM(xyz_pairs_first),get_COM(xyz_pairs_last)) ]
rot = rotation_matrix_align_to_x_axis(vec1)

print('Rotation : ' + str(rot))

# write pdb files
f_chain = PdbFile(outfilename)
f_chain.open_to_write()

for c in chains :
    for r in c.residues :
        for a in r.atoms :
            a.xyz.transform_rot(rot) 

f_chain.write_all(chains)
f_chain.close()     


print('***ALL DONE NORMALLY***')




