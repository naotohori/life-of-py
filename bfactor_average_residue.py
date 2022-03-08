#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
@author: Naoto Hori
'''
import sys
from lop.file_io.pdb import PdbFile

if len(sys.argv) != 4:
    print ('\n Usage: SCRIPT [bfactor file] [PDB file] [output bfactor file]\n')
    sys.exit(2)

f_bf_in = open(sys.argv[1], 'r')
f_pdb = PdbFile(sys.argv[2])
f_pdb.open_to_read()
chains = f_pdb.read_all()
f_pdb.close()
f_bf_out = open(sys.argv[3], 'w')

f_bf_out.write('## SCRIPT: bfactor_average_residue.py\n')
f_bf_out.write('## argv[1]: '+sys.argv[1]+'\n')
f_bf_out.write('## argv[2]: '+sys.argv[2]+'\n')
f_bf_out.write('## argv[3]: '+sys.argv[3]+'\n')
f_bf_out.write('\n\n')

bf = []
for line in f_bf_in :
    if line.find('#') != -1 :
        continue
    bf.append(float(line.split()[1]))
    
# current version is only for single chain
#if len(chains) != 1:
#    print ('Error: len(chains) != 1')
#    sys.exit(2)
    
imp = 0
for c in chains:
    sum_np = 0.0
    n_np = 0
    for r in c.residues:
        sum_p = sum_np
        sum_b = 0.0
        sum_s = 0.0
        sum_np = 0.0
        sum_pro = 0.0
        n_p = n_np
        n_b = 0
        n_s = 0
        n_np = 0
        n_pro = 0
        for a in r.atoms:
            if a.name[0:1] == 'H' or a.name[0:2] == ' H' :
    #            print a.name, 'Hydrogen'
                pass
            elif a.name == ' CA ':
                sum_pro += bf[imp]
                n_pro += 1
            elif a.name == " O3'" :
    #            print a.name, 'next P'
                sum_np += bf[imp] 
                n_np += 1
            elif a.name.find("'") != -1:
    #            print a.name, 'S'
                sum_s += bf[imp]
                n_s += 1
            elif a.name in (" P  ", " OP1", " OP2", "OP3") :
    #            print a.name, 'P'
                sum_p += bf[imp]
                n_p += 1
            else :
    #            print imp,a.name, 'B'
                sum_b += bf[imp]
                n_b += 1
            imp += 1
        if n_pro == 0:
            if n_p != 0:
                f_bf_out.write('%f\n' % (sum_p/float(n_p) ,))
            if n_s != 0:
                f_bf_out.write('%f\n' % (sum_s/float(n_s) ,))
            if n_b != 0:
                f_bf_out.write('%f\n' % (sum_b/float(n_b) ,))
        else:
            f_bf_out.write('%f\n' % (sum_pro/float(n_pro), ))
