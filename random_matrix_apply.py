#!/usr/bin/env python

import sys
import os
import glob
import re

if len(sys.argv) != 4:
    print('Usage: % SCRIPT [pdb file] [matrix DIR] [output prefix]')
    sys.exit(2)

filename_pdb = sys.argv[1]
dir_matrix = sys.argv[2]
prefix_out = sys.argv[3]

#matfiles = glob.glob(dir_matrix+'*')
#nfile = len(matfiles)

for file in glob.glob(dir_matrix+'*') :
     id = re.match(".*\/_(\d*).mat",file).group(1)
     filename_new = prefix_out + id + '.pdb'
     os.system("bestfit_pdb_by_matrix.py "+filename_pdb+" "+file+" "+filename_new)

