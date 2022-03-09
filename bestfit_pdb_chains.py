#!/usr/bin/env python

from lop.file_io.pdb import PdbFile
from CalcROT import calcrotation
from lop.util_pdb import chains_to_ndarray
from numpy import dot

def fit(chains_ref, chains_que):
    d_ref = chains_to_ndarray(chains_ref)
    d_que = chains_to_ndarray(chains_que)

    rmsd, rot = calcrotation(d_ref.T, d_que.T)

    return rmsd, rot 

def apply_rot(chains, rot):

    for c in chains:
        for r in c.residues:
            for a in r.atoms:
                coords = a.xyz.get_as_list()
                coords = dot(rot, coords+[1.0,])[0:3]
                a.xyz.x = coords[0]
                a.xyz.y = coords[1]
                a.xyz.z = coords[2]


if __name__ == '__main__':

    import sys
    import argparse

    parser = argparse.ArgumentParser(
             description='Calculate RMSD and transformation matrix',
             formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('pdb_ref', help='input reference PDB')
    parser.add_argument('pdb_que', help='input query PDB')
    parser.add_argument('--mtx', dest='file_mtx', help='output matrix (mtx) file')
    args = parser.parse_args()

    chains_ref = PdbFile(sys.argv[1],'r').read_all_and_close()
    chains_que = PdbFile(sys.argv[2],'r').read_all_and_close()

    rmsd, rot = fit(chains_ref, chains_que)

    print(rmsd)
    print(rot)

    if args.file_mtx is not None:
        f_out = open(args.file_mtx,'w')
        f_out.write('#matrix\n')
        for i in range(4):
            for j in range(4):
                f_out.write(' %f' % (rot[i,j],))
            f_out.write('\n')

