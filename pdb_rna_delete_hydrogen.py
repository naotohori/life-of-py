#!/usr/bin/env python

import os
import sys

sys.argv
if len(sys.argv) != 3:
    print() 
    print('This script deletes hydrogen atoms from PDB file.')
    print('Usage: % SCRIPT [PDB file] [output PDB file]')
    print() 
    sys.exit(2)

filename_pdb = sys.argv[1]
filename_out = sys.argv[2]
file_out = file(filename_out,'w')

file_out.write('REMARK script: '+sys.argv[0]+'\n')
file_out.write('REMARK cwd   : '+os.getcwd()+'\n')
file_out.write('REMARK input : '+filename_pdb+'\n')
file_out.write('REMARK output: '+filename_out+'\n')

for line in open(filename_pdb) :

    if line[0:4] == 'ATOM' or line[0:6] == 'HETATM' :
        if line[12:13] != 'H' and line[12:14] != ' H' :
            file_out.write(line)

file_out.write('TER\n')
file_out.close()

