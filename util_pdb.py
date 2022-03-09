'''
@author: Naoto Hori
'''

from numpy import zeros, float64

def chains_to_ndarray(chains):

    natom = 0
    for c in chains:
        natom += c.num_atom()

    d = zeros((natom, 3), dtype=float64, order='C')
    iatom = 0
    for c in chains:
        for r in c.residues:
            for a in r.atoms:
                d[iatom] = a.xyz.get_as_ndarray()
                iatom += 1
    return  d
