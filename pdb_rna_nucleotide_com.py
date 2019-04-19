#!/usr/bin/env python

'''
Read all-atom pdb
Return pdb which contains center of masses of nucleotides.
(ignore Hydrogen atom)
'''

import sys
from cafysis.file_io.pdb import PdbFile
from cafysis.elements.coord import Coord
from cafysis.elements.pdb import Chain, Residue, Atom


element2mass = {'P':30.973761, 'O':15.9994, 'C':12.0107, 'N':14.0067}

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: %SCRIPT [input PDB file] [output PDB]')
        sys.exit(2)
        

pdb_in = PdbFile(sys.argv[1])
pdb_in.open_to_read()
chains = pdb_in.read_all()
pdb_in.close()

pdb_out = PdbFile(sys.argv[2])
pdb_out.open_to_write()
pdb_out.write_remark('Generated using cafysis/pdb_rna_nucleotide_com.py')
import time
import datetime
pdb_out.write_remark('At %s' % (datetime.datetime.fromtimestamp(time.time()).strftime('%Y%m%d %H:%M%S'),))

chains_com = []

from . import elements

i_serial = 0

for c in chains:
    c_com = Chain()
    for r in c.residues:
        a_com = Atom()
        i_serial += 1
        a_com.serial = i_serial
        a_com.name = 'NUC'
        a_com.res_name = r.atoms[0].res_name
        a_com.chain_id = r.atoms[0].chain_id
        a_com.res_seq = r.atoms[0].res_seq

        mass_sum = 0.0
        xyz_sum = Coord()
        for a in r.atoms:
            if a.element.strip() == 'H':
                continue
            if a.element.strip() in element2mass:
                mass = element2mass[a.element.strip()]
            else:
                print('Error: no key %s in element2mass.' % (a.element,))
            xyz_sum += a.xyz * mass
            mass_sum += mass
        a_com.xyz = xyz_sum / float(mass_sum)
        
        r_com = Residue()
        r_com.push_atom(a_com)
        c_com.push_residue(r_com)
    
    chains_com.append(c_com)

pdb_out.write_all(chains_com)

