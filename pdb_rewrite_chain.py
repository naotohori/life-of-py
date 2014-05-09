#!/usr/bin/env python

import os
import sys

if len(sys.argv) != 4:
    print 
    print 'Usage: % SCRIPT [PDB file] [char]  [output PDB file]'
    print 
    sys.exit(2)

filename_pdb = sys.argv[1]
chain = sys.argv[2]
filename_out = sys.argv[3]
file_out = file(filename_out,'w')

file_out.write('REMARK script: '+sys.argv[0]+'\n')
file_out.write('REMARK cwd   : '+os.getcwd()+'\n')
file_out.write('REMARK input : '+filename_pdb+'\n')
file_out.write('REMARK char  : '+chain+'\n')
file_out.write('REMARK output: '+filename_out+'\n')

for line in open(filename_pdb) :

    if line[0:4] == 'ATOM' or line[0:6] == 'HETATM' :
            new_line = line[:21] + chain + line[22:]
            file_out.write(new_line)

file_out.write('TER\n')
file_out.close()

