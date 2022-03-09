#!/usr/bin/env python

from lop.elements.pdb import Chain, Residue, Atom
from lop.elements.coord import Coord
from lop.para.mass import ATOM_MASS

import sys

## If residues in the original PDB are
# DT5, DA, DC ...., then
SEQ_POSITION = 1
# T5, A, C ...., then
#SEQ_POSITION = 0

ATOMS_P = ('P', 'OP1', 'OP2', "O5'")  # "O3'" from the previous nucleotide
ATOMS_P_AMBER = ("P", "O1P", "O2P", "O5*")
ATOMS_S = ("C5'", "C4'", "C3'", "C2'", "C1'", "O4'")
ATOMS_S_AMBER = ("C5*", "C4*", "C3*", "C2*", "C1*", "O4*")
ATOMS_A = ("N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4")
ATOMS_G = ("N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4")
ATOMS_C = ("N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6")
ATOMS_T = ("N1", "C2", "O2", "N3", "C4", "O4", "C5", "C7", "C6")


def aa2cg_dna(chains, flg_start_P=False):

    cg_chains = []

    atom_id = 0
    for c in chains:
        c_cg = Chain()
        xyz_O3 = None
        res_id = 0

        for ir, r in enumerate(c.residues):

            xyz_P = Coord()
            nP = 0

            if xyz_O3 is None:
                if ir != 0:
                    print("Error: cannot find O3 from the previous nucleotide")
                    sys.exit(2)
            else:
                xyz_P += xyz_O3 * ATOM_MASS["O"]
                nP += ATOM_MASS["O"]
                xyz_O3 = None

            xyz_S = Coord()
            nS = 0 
            xyz_B = Coord()
            nB = 0

            for a in r.atoms:
                name = a.name.strip()

                if name[0] == 'H':
                    continue

                elif (name == "O3'") or (name == "O3*"):
                    xyz_O3 = a.xyz

                elif (name in ATOMS_P) or (name in ATOMS_P_AMBER) :
                    if name.find("O") != -1:
                        xyz_P += a.xyz * ATOM_MASS["O"]
                        nP += ATOM_MASS["O"]
                    elif name[0] == 'P':
                        xyz_P += a.xyz * ATOM_MASS["P"]
                        nP += ATOM_MASS["P"]

                elif (name in ATOMS_S) or (name in ATOMS_S_AMBER) :
                    if name.find("C") != -1:
                        xyz_S += a.xyz * ATOM_MASS["C"]
                        nS += ATOM_MASS["C"]
                    else:
                        xyz_S += a.xyz * ATOM_MASS["O"]
                        nS += ATOM_MASS["O"]

                else:
                    if name.find("N") != -1:
                        xyz_B += a.xyz * ATOM_MASS["N"]
                        nB += ATOM_MASS["N"]
                    elif name.find("C") != -1:
                        xyz_B += a.xyz * ATOM_MASS["C"]
                        nB += ATOM_MASS["C"]
                    elif name.find("O") != -1:
                        xyz_B += a.xyz * ATOM_MASS["O"]
                        nB += ATOM_MASS["O"]

                nt = a.res_name.strip()
                nchain = a.chain_id.strip()

            res_id += 1
            r_cg = Residue()

            if ir != 0:
                atom_id += 1
                a = Atom()
                a.serial = atom_id
                a.name = ' P  '
                a.res_name = 'D%s ' % nt[SEQ_POSITION:SEQ_POSITION+1]
                a.chain_id = '%s' % nchain
                a.res_seq = res_id
                a.xyz = xyz_P / float(nP)
                r_cg.push_atom(a)

            atom_id += 1
            a = Atom()
            a.serial = atom_id
            a.name = ' S  '
            a.res_name = 'D%s ' % nt[SEQ_POSITION:SEQ_POSITION+1]
            #a.chain_id = 'A'
            a.chain_id = '%s' % nchain
            a.res_seq = res_id
            a.xyz = xyz_S / float(nS)
            r_cg.push_atom(a)

            atom_id += 1
            a = Atom()
            a.serial = atom_id
            #a.name = ' %sb '  %  nt
            a.name = ' %sb '  %  nt[SEQ_POSITION:SEQ_POSITION+1]
            #a.res_name = 'R%s ' % nt
            a.res_name = 'D%s ' % nt[SEQ_POSITION:SEQ_POSITION+1]
            #a.chain_id = 'A'
            a.chain_id = '%s' % nchain
            a.res_seq = res_id
            a.xyz = xyz_B / float(nB)
            r_cg.push_atom(a)

            c_cg.push_residue(r_cg)

        cg_chains.append(c_cg)

    return cg_chains


if __name__ == '__main__':

    from cafysis.file_io.pdb import PdbFile

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

    chains_cg = aa2cg_dna(chains, flg_start_P)

    cg = PdbFile(sys.argv[-1], 'w')
    cg.write_all(chains_cg)
    cg.close()

