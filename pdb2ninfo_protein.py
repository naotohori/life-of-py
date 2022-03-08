#!/usr/bin/env python

'''
Created on 2015/06/30
@author: Naoto Hori
'''

import sys
from lop.elements.ninfo import Fene, LJ

fene_SC = False #  True: 2-beads model,  False: 1-bead model
fene_coef  = 20.0
fene_dist2 =  2.0

LJ_cutoff = 8.0
LJ_min_ij = 3  # j-i >= 3
LJ_coef = 2.0

def generate_FENEs(chains, ns):
    iunit = 0
    imp_offset = 0
    for c in chains:
        iunit += 1
        for i in range(c.num_res()-1):
            if fene_SC:
                print ('Error: side chain modeis has not been implemented yet!!')
                sys.exit(2)
            else:
                ninfo = Fene(iunit1=iunit, iunit2=iunit, imp1=imp_offset+i+1, imp2=imp_offset+i+2, imp1un=i+1, imp2un=i+2)
                xyzi = c.residues[i].find_Calpha_atom().xyz
                ninfo.native = c.residues[i+1].find_Calpha_atom().xyz.distance( xyzi )
                ninfo.coef = fene_coef
                ninfo.dist2 = fene_dist2
                ns.fenes.append(ninfo)
        imp_offset += c.num_res()
    
    return


def generate_LJs(chains, ns):

    imp1_offset = 0
    for ic1, c1 in enumerate(chains):

        # intra chain
        for ir1 in range(c1.num_res()):
            xyz1 = c1.residues[ir1].find_Calpha_atom().xyz

            for ir2 in range(ir1+LJ_min_ij, c1.num_res()):
                d = c1.residues[ir2].find_Calpha_atom().xyz.distance( xyz1 )
                if d <= LJ_cutoff:
                    ninfo = LJ(iunit1 = ic1+1, iunit2 = ic1+1, 
                               imp1 = imp1_offset+ir1+1, imp2 = imp1_offset+ir2+1, 
                               imp1un = ir1+1, imp2un = ir2+1)
                    ninfo.native = d
                    ninfo.coef = LJ_coef
                    ns.LJs.append(ninfo)

        # inter chain
        imp2_offset = 0
        for ic2 in range(ic1+1, len(chains)):
            c2 = chains[ic2]

            for ir1 in range(c1.num_res()):
                xyz1 = c1.residues[ir1].find_Calpha_atom().xyz

                for ir2 in range(c2.num_res()):
                    d = c2.residues[ir2].find_Calpha_atom().xyz.distance( xyz1 )
                    if d<= LJ_cutoff:
                        ninfo = LJ(iunit1 = ic1+1, iunit2 = ic2+1, 
                               imp1 = imp1_offset+ir1+1, imp2 = imp2_offset+ir2+1, 
                               imp1un = ir1+1, imp2un = ir2+1)
                        ninfo.native = d
                        ninfo.coef = LJ_coef
                        ns.LJs.append(ninfo)
            imp2_offset += c2.num_res()
            
        imp1_offset += c1.num_res()
    return


if __name__ == '__main__':
    import sys
    if len(sys.argv) != 3:
        print ('Usage: %SCRIPT [input PDB] [output ninfo]')
        sys.exit(2)

    from lop.file_io.pdb import PdbFile
    from lop.file_io.ninfo import NinfoFile
    from lop.elements.ninfo import NinfoSet

    f = PdbFile(sys.argv[1])
    f.open_to_read()
    chains = f.read_all()
    f.close()

    ns = NinfoSet()
    generate_FENEs(chains, ns)
    generate_LJs(chains, ns)

    f_ninfo = NinfoFile(sys.argv[-1])
    f_ninfo.open_to_write()
    f_ninfo.write_all(ns)
    f_ninfo.close()

