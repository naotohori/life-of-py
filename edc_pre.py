#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 2011/05/29 coded by Naoto HORI

from .file_io.pdb import PdbFile
from .file_io.ninfo import NinfoFile
from .elements.ninfo import NinfoSet
import sys

if len(sys.argv) < 9 or len(sys.argv)%2 != 1 :
    print(('Usage: SCRIPT [input PDB] [cafe_to_amber file] [ninfo file] '+
           '[(imp begin, imp end) ...] [log file] [Amber file] [residue file]'))
    sys.exit(2)
    
# Output file
f_out     = open(sys.argv[-3], 'w')
f_amber   = open(sys.argv[-2], 'w')
f_residue = open(sys.argv[-1], 'w')
for i,arg in enumerate(sys.argv):
    f_out.write('argv[%i] = %s\n' %  (i+1, arg))
f_out.write('\n')
    
# Reading cafe_to_amber
aa_res = {}  # put in dictionary
for line in open(sys.argv[2]) :
    if line.find('#') != -1 :
        continue
    linesp = line.split()
    # Read data
    ichain   = int(linesp[0].strip())   
    ires_ini = int(linesp[3].strip())
    ires_end = int(linesp[4].strip())
    aa_res[ichain] = (ires_ini, ires_end)

# Reading cafemol PDB
pdb_cafemol = PdbFile(sys.argv[1])
pdb_cafemol.open_to_read()
cafe_chains = pdb_cafemol.read_all()
pdb_cafemol.close()

# Generate "cafe_mp_info" data
offset_res = 0
imp = 0
cafe_mp_info = {}
for (ic, c) in enumerate(cafe_chains):
    for (ir, r) in enumerate(c.residues):
        for (ia, a) in enumerate(r.atoms):
            imp += 1
            cafe_mp_info[imp] = (ic+1,   # chain
                                 a.res_seq+offset_res,  # res
                                 a.res_seq,  # res(local)
                                 ia+1,       # imp(local)
                                 a.res_name,
                                 a.name)
            f_out.write('%6i %3i %6i %6i %s %6i %s\n'
                         % (imp, ic+1, a.res_seq+offset_res, a.res_seq,
                            a.res_name, a.res_seq, a.name))
    offset_res += c.get_atom(-1).res_seq
f_out.write('\n')

nchain = len(cafe_chains)
f_out.write('#nchain = %i\n\n' % nchain)

# Generate "cafe_res"
cafe_res = {}
sum_res = 0
for (i,c) in enumerate(cafe_chains):
    nres = c.get_atom(-1).res_seq
    cafe_res[i+1] = (sum_res+1, sum_res+nres)
    sum_res += nres

# check
f_out.write('#check\n')
for i in range(len(cafe_res)) :
    n_amber = aa_res[i+1][1] - aa_res[i+1][0] + 1
    n_cafe = cafe_res[i+1][1]- cafe_res[i+1][0] + 1
    if n_cafe != n_amber :
        print('Inconsistent!!')
        print(i+1, 'amber=', n_amber, 'cafe=', n_cafe)
    f_out.write('%3i , amber= %6i-%6i , cafe= %6i-%6i\n' % (i+1, aa_res[i+1][0], aa_res[i+1][1],
                                        cafe_res[i+1][0], cafe_res[i+1][1]) )
f_out.write('\n\n')

# Reading .ninfo File
f_ninfo = NinfoFile(sys.argv[3])
f_ninfo.open_to_read()
ninfo = NinfoSet()
f_ninfo.read_all(ninfo)

# Reading arguments (ID begin, ID end)...
id_pairs =[]
for (i,arg) in enumerate(sys.argv[4:-3]) :
    if i%2 == 0:
        pair = (int(arg), )
    else :
        id_pairs.append( pair + (int(arg), ))
#print id_pairs

# Generate ligand set : imp(cafe)
ligand_set = set()
for pair in id_pairs :
    for imp in range(pair[0], pair[1]+1) :
        ligand_set.add(imp)
f_out.write('#Ligands imp(cafe)\n')
for i in ligand_set :
    f_out.write('%i ' % i)
f_out.write('\n\n')

# Generate recepter set by searching contacts of ninfo : imp(cafe)
#  (contact and basepair)
recepter_set = set()
for lmp in ligand_set :
    for con in ninfo.contacts :
        if lmp == con.imp1 :
            recepter_set.add(con.imp2)
        elif lmp == con.imp2 :
            recepter_set.add(con.imp1)
    for bp in ninfo.basepairs :
        if lmp == bp.imp1 :
            recepter_set.add(bp.imp2)
        elif lmp == bp.imp2 :
            recepter_set.add(bp.imp1)
f_out.write('#Recepters imp(cafe) BEFORE removing dumplication\n')
for i in recepter_set :
    f_out.write('%i ' % i)
f_out.write('\n\n')
            
# Remove duplicate imp from recepters : imp(cafe)
for lmp in ligand_set :
    recepter_set.discard(lmp)
f_out.write('#Recepters imp(cafe) AFTER removing dumplication\n')
for i in recepter_set :
    f_out.write('%i ' % i)
f_out.write('\n\n')

# convert from cafemol-imp to amber-ires
## ligand : ires(amber)
ligand_output_pairs = []
for pair in id_pairs:
    ichain = cafe_mp_info[pair[0]][0] # chain
    ires = cafe_mp_info[pair[0]][1] # ires(cafe)
    #chain = 0 
    #for ichain in xrange(1, nchain+1) :
    #    if ires >= cafe_res[ichain][0] and ires <= cafe_res[ichain][1] :
    #        chain = ichain
    #        break
    #begin =  aa_res[chain][0] + ires - cafe_res[chain][0]  # ires(amber)
    begin =  aa_res[ichain][0] + ires - cafe_res[ichain][0]  # ires(amber)
    
    ichain = cafe_mp_info[pair[1]][0]
    ires = cafe_mp_info[pair[1]][1]
    #chain = 0 
    #for ichain in xrange(1, nchain+1) :
    #    if ires >= cafe_res[ichain][0] and ires <= cafe_res[ichain][1] :
    #        chain = ichain
    #        break
    #end =  aa_res[chain][0] + ires - cafe_res[chain][0]   # ires(amber)
    end =  aa_res[ichain][0] + ires - cafe_res[ichain][0]   # ires(amber)
    ligand_output_pairs.append( (begin, end) )

## recepter
recepter_set_amber = set()
for imp in recepter_set :
    ires = cafe_mp_info[imp][1]  # ires(cafe)
    ichain = cafe_mp_info[imp][0]  # chain
    #chain = 0
    #for ichain in xrange(1, nchain+1) :
    #    if ires >= cafe_res[ichain][0] and ires <= cafe_res[ichain][1] :
    #        chain = ichain
    #        break
    #recepter_set_amber.add( aa_res[chain][0] + ires - cafe_res[chain][0] )
    recepter_set_amber.add( aa_res[ichain][0] + ires - cafe_res[ichain][0] )
    
f_out.write('#Recepters ires(amber)\n')
for i in recepter_set_amber :
    f_out.write('%i ' % i)
f_out.write('\n\n')

recepters = list(recepter_set_amber)
recepters.sort(cmp=None, key=None, reverse=False)
f_out.write('#Recepters ires(amber) sorted\n')
for i in recepters :
    f_out.write('%i ' % i)
f_out.write('\n\n')

recepter_output_pairs = []
imp_begin = recepters[0]
imp_pre = imp_begin
for imp in recepters[1:] :
    if imp != imp_pre + 1 :
        recepter_output_pairs.append( (imp_begin, imp_pre ) )
        imp_begin = imp
        imp_pre = imp
    else:
        imp_pre = imp
#if recepters[-2] + 1 == recepters[-1] :
recepter_output_pairs.append( (imp_begin, imp_pre) )
    
### Residue file output
f_residue.write('#ichain ires  imp ires_l  imp_l  res  mp   ires_aa\n')
for imp in range(1, len(cafe_mp_info)+1) :
    mp_info = cafe_mp_info[imp]
    ichain  = cafe_mp_info[imp][0]
    ires    = cafe_mp_info[imp][1]
    ires_l  = cafe_mp_info[imp][2]
    imp_l   = cafe_mp_info[imp][3]
    res_name= cafe_mp_info[imp][4]
    mp_name = cafe_mp_info[imp][5]
    ires_amber = aa_res[ichain][0] + ires - cafe_res[ichain][0] 
    f_residue.write('%3i %6i %6i %6i %6i %4s %4s %6i\n' 
                % (ichain, ires, imp, ires_l, imp_l, res_name, mp_name, ires_amber))
    
### Amber file output
#f_amber.write('############## OUTPUT ################\n')
f_amber.write('Ligands\n')
for pair in ligand_output_pairs :
    f_amber.write('LRES %i %i\n' % pair)
f_amber.write('END\n')
f_amber.write('Recepters\n')
for pair in recepter_output_pairs :
    f_amber.write('RRES %i %i\n' % pair)
f_amber.write('END\n')
f_amber.write('Residues to print\n')
for pair in ligand_output_pairs :
    f_amber.write('RES %i %i\n' % pair)
for pair in recepter_output_pairs :
    f_amber.write('RES %i %i\n' % pair)
f_amber.write('END\n')
f_amber.write('END\n')
            