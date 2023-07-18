#!/usr/bin/env python
# Huong May 10th 2023, Guangzhou

from lop.elements.pdb import Chain, Residue, Atom
from lop.elements.coord import Coord

import sys

def change_name_to_C4(chains):
    chains_cg = []
    for ic, c in enumerate(chains):
        c_cg = Chain() 
        for r in c.residues:
            r_cg = Residue()

            for a in r.atoms:
                a.name = "C4\'"
                r_cg.push_atom(a)

            c_cg.push_residue(r_cg)
    
        chains_cg.append(c_cg)

    return chains_cg

def change_res_name(chains):
    chains_cg = []
    for ic, c in enumerate(chains):
        c_cg = Chain() 
        for r in c.residues:
            r_cg = Residue()

            for a in r.atoms:
                a.res_name = a.name.strip()
                r_cg.push_atom(a)

            c_cg.push_residue(r_cg)
    
        chains_cg.append(c_cg)

    return chains_cg

def change_resid_to_fasta(chains, fasta_file):

    flg_new_chain = False
    seq_list = []
    seq = ''
    for l in open(fasta_file,'r'):
        if l[0] == '>' :
            if flg_new_chain == True:
                seq_list.append(seq)
            flg_new_chain = True
            seq = ''
        elif l[0] == '#' or l[0] == ';':
            continue
        else:
            seq += l.strip()
    
    seq_list.append(seq)
    
    print(seq_list)
    
#    if len(flatten(seq_list)) != len(chains):
#        print('Number of chains in fasta file '+str(len(seq_list))+' does not match with number of chains in pdb file ' +str(len(chains)))
#    n_nt = len(flatten(seq_list))

    chains_cg = []
    ic = -1
    ir = -1
    chain_id = ''
    res_name
    c_cg = Chain()
    for c in chains:
        
        for r in c.residues:
        
            r_cg = Residue()
            if chain_id!=r.atoms[0].chain_id:
                if chain_id!='':
                    if ir!= len(seq_list[ic]):
                        print('Warning: number of residues '+str(len(seq_list[ic]))+'in chain '+str(ic)+' in fasta file does not match '+str(ir)+' got from pdb file')
                    chains_cg.append(c_cg)
                    c_cg = Chain() 
                ic +=1
                chain_id=r.atoms[0].chain_id
                
            if r.atoms[0].res_name.strip() != seq_list[ic][ir]:
                print('Chain '+ str(ic) +' residue '+ str(ir) +' sequence '+ r.atoms[0].res_name.strip() +' does not match with fasta sequence '+seq_list[ic][ir])
                return(2)

            for a in r.atoms:
                ir+=1
                a.res_seq = ir+1
                r_cg.push_atom(a)

            c_cg.push_residue(r_cg)
    
    
    chains_cg.append(c_cg)

    return chains_cg

def change_resid_to_index(chains):
    chains_cg = []
    for ic, c in enumerate(chains):
        c_cg = Chain() 
        for r in c.residues:
            r_cg = Residue()

            for a in r.atoms:
                a.res_seq = a.serial
                print(a.res_seq)
                r_cg.push_atom(a)

            c_cg.push_residue(r_cg)
    
        chains_cg.append(c_cg)

    return chains_cg
    


if __name__ == '__main__':
    
    from lop.file_io.pdb import PdbFile

    if len(sys.argv) not in (3,4):
        print('Usage: SCRIPT [input cg PDB VMD wrote from xyz] [fasta file (optional)] [output cg PDB]')
        sys.exit(2)

    aa = PdbFile(sys.argv[1], 'r')
    chains = aa.read_all()
    aa.close()

    chains_cg = change_res_name(chains)
    chains_cg = change_name_to_C4(chains_cg)

    if len(sys.argv) == 4:
        fasta_file = sys.argv[2]

        chains_cg = change_resid_to_fasta(chains_cg, fasta_file)
    else:
        chains_cg = change_resid_to_index(chains_cg)


    cg = PdbFile(sys.argv[-1], 'w')
    cg.write_all(chains_cg)
    cg.close()


