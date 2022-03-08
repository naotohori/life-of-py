#!/usr/bin/env python

from lop.elements.pdb import Chain, Residue, Atom
from lop.elements.coord import Coord

import sys

ATOMS_P = ('P', 'OP1', 'OP2')

def aa2cg(chains_aa, flg_start_P=False):

    chains_cg = []
    
    res_id = 0
    atom_id = 0
    for c in chains_aa:
        c_cg = Chain()
        #xyz_O3 = None
        flg_noS = False
        flg_noB = False

        for ir, r in enumerate(c.residues):

            # Check if the previous residue was normal
            if flg_noB:
                print("Error: there was no atoms for B, ir =", ir-1)
                sys.exit(2)
            if flg_noS:
                print("Error: there was no atoms for S, ir =", ir-1)
                sys.exit(2)

            xyz_P = Coord()
            nP = 0
            #if not xyz_O3 is None:
            #    xyz_P += xyz_O3
            #    nP += 1
            #    xyz_O3 = None
            xyz_S = Coord()
            nS = 0 
            xyz_B = Coord()
            nB = 0
    
            for a in r.atoms:
                name = a.name.strip()
                if name[0] == 'H':
                    continue
                #elif name == "O3'":
                #    xyz_O3 = a.xyz
                elif name in ATOMS_P:
                    xyz_P += a.xyz
                    nP += 1
                elif name.find("'") != -1:
                    xyz_S += a.xyz
                    nS += 1
                else:
                    xyz_B += a.xyz
                    nB += 1
                nt = a.res_name.strip()
    
            res_id += 1
            r_cg = Residue()
    
            if ir != 0 or flg_start_P:
                if nP == 0:
                    print('Error: no atoms for P, ir =', ir)
                    sys.exit(2)
                atom_id += 1
                a = Atom()
                a.serial = atom_id
                a.name = ' P  '
                a.res_name = 'R%s ' % nt
                a.chain_id = 'A'
                a.res_seq = res_id
                a.xyz = xyz_P / float(nP)
                r_cg.push_atom(a)
    
            if nS > 0:
                atom_id += 1
                a = Atom()
                a.serial = atom_id
                a.name = ' S  '
                a.res_name = 'R%s ' % nt
                a.chain_id = 'A'
                a.res_seq = res_id
                a.xyz = xyz_S / float(nS)
                r_cg.push_atom(a)
            else:
                flg_noS = True
    
            if nB > 0:
                atom_id += 1
                a = Atom()
                a.serial = atom_id
                a.name = ' %sb '  %  nt
                a.res_name = 'R%s ' % nt
                a.chain_id = 'A'
                a.res_seq = res_id
                a.xyz = xyz_B / float(nB)
                r_cg.push_atom(a)
            else:
                flg_noB = True
    
            c_cg.push_residue(r_cg)
    
        chains_cg.append(c_cg)

    return chains_cg


if __name__ == '__main__':
    
    from lop.file_io.pdb import PdbFile

    if len(sys.argv) not in (3,4):
        print('Usage: SCRIPT [input aa PDB] [flag to start with P (Y/n)] [output cg PDB]')
        sys.exit(2)

    flg_start_P = False
    if len(sys.argv) == 4:
        if sys.argv[2] == 'Y':
            flg_start_P = True


    aa = PdbFile(sys.argv[1], 'r')
    chains = aa.read_all()
    aa.close()

    chains_cg = aa2cg(chains, flg_start_P)

    cg = PdbFile(sys.argv[-1], 'w')
    cg.write_all(chains_cg)
    cg.close()

