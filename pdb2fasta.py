#!/usr/bin/env python

def show_usage() :
    print('')
    print(' usage: SCRIPT [input pdb] [output fasta] [[TITLE]]')
    print('')
    
import sys
if not len(sys.argv) in (3,4) :
    show_usage()
    sys.exit(2)

from cafysis.file_io.pdb import PdbFile

f_pdb = PdbFile(sys.argv[1])
f_pdb.open_to_read()
chains = f_pdb.read_all()

f_fasta = open(sys.argv[2],'w')

if len(sys.argv) == 4 :
    f_fasta.write('>%s|PDBID|CHAIN|SEQUENCE' % sys.argv[3] + '\n')
else :
    f_fasta.write('>%s|PDBID|CHAIN|SEQUENCE' % sys.argv[1][0:6] + '\n')
    
amino = {'ALA':'A', 'ASX':'B', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F',
         'GLY':'G', 'HIS':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L', 'MET':'M',
         'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R', 'SER':'S', 'THR':'T',
         'SEC':'U', 'VAL':'V', 'TRP':'W', 'XAA':'X', 'TYR':'Y', 'GLX':'Z', }
rna = { '  A':'A', '  U':'U', '  G':'G', '  C':'C',
        ' RA':'A', ' RU':'U', ' RG':'G', ' RC':'C',
        'RA ':'A', 'RU ':'U', 'RG ':'G', 'RC ':'C' }

str = ''
for c in chains :
    for r in c.residues :
        c3 = r.atoms[0].res_name
        if c3 in amino:
            c1 = amino[c3]
        elif c3 in rna :
            c1 = rna[c3]
        else:
            print('no data', c3)
            sys.exit(2)
        str += c1
        
#print str

lines = int (len(str) / 80)
for i in range(lines) :
    f_fasta.write(str[i*80:(i+1)*80] + '\n')
f_fasta.write(str[lines*80:] + '\n')
