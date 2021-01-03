#!/usr/bin/env python

from cafysis.file_io.pdb import PdbFile
from CalcROT import calcrotation
from cafysis.util_pdb import chains_to_ndarray
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
    if len(sys.argv) != 3:
        print ("Usage: SCRIPT [input reference PDB] [input query PDB]")
        sys.exit(2)

    chains_ref = PdbFile(sys.argv[1],'r').read_all_and_close()
    chains_que = PdbFile(sys.argv[2],'r').read_all_and_close()

    rmsd, rot = fit(chains_ref, chains_que)

    print(rmsd)
    print(rot)

