#!/usr/bin/env python

import os
import sys

sys.argv
if len(sys.argv) != 3:
    print 
    print 'This script extracts backbone-atoms which have name "?OP?","?P??","???\'" from PDB file.'
    print 'Usage: % SCRIPT [PDB file] [output PDB file]'
    print 
    sys.exit(2)

filename_pdb = sys.argv[1]
filename_out = sys.argv[2]
file_out = file(filename_out,'w')

file_out.write('REMARK script: '+sys.argv[0]+'\n')
file_out.write('REMARK cwd   : '+os.getcwd()+'\n')
file_out.write('REMARK input : '+filename_pdb+'\n')
file_out.write('REMARK output: '+filename_out+'\n')

pre_res_id = 0
backbone = [None, None, None, None, None, None]
for line in open(filename_pdb) :

    if line[0:4] == 'ATOM' or line[0:6] == 'HETATM' :
        res_id = int(line[22:26])
        if pre_res_id == 0:
            pre_res_id = res_id

        if res_id != pre_res_id :
            for i in xrange(6) :
                if backbone[i] :
                    file_out.write(backbone[i])
            backbone = [None, None, None, None, None, None]
            pre_res_id = res_id

        if line[12:16] == " P  " : 
            backbone[0] = line
        if line[12:16] == " O5'" :
            backbone[1] = line
        if line[12:16] == " C5'" :
            backbone[2] = line
        if line[12:16] == " C4'" :
            backbone[3] = line
        if line[12:16] == " C3'" :
            backbone[4] = line
        if line[12:16] == " O3'" :
            backbone[5] = line

for i in xrange(6) :
    if backbone[i] : 
        file_out.write(backbone[i])

file_out.write('TER\n')
file_out.close()

