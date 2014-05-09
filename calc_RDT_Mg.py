#!/usr/bin/env python
# -*- coding: utf-8 -*-

from file_pdb import PdbFile
import sys
from numpy import histogram

if len(sys.argv) < 6:
    print ('\n Usage: % SCRIPT [min] [max] [width] [PDB] [[PDB] ....] [output]\n')
    sys.exit(2)
    
# parameters
HIST_MIN = float(sys.argv[1])
HIST_MAX = float(sys.argv[2])
HIST_WID = float(sys.argv[3])
    
# Read chains from all input files
chains = []
for arg in sys.argv[4:-1] :
    file = PdbFile(arg)
    file.open_to_read()
    file.flg_HETATM = True
    chains += file.read_all()
    
# Output file preparing
import os
if os.path.exists(sys.argv[-1]) :
    if raw_input('File %s exists. Rewrite? [Y/n]'%sys.argv[-1]) != 'Y' :
        sys.exit(2)
file_out = open(sys.argv[-1], 'w')

# Search Mg
MgIds = []
for (i_chain,chain) in enumerate(chains) :
    for (i_res, res) in enumerate(chain.residues) :
        for (i_atom, atom) in enumerate(res.atoms) :
            #print (i_chain,i_res,i_atom,atom.name)
            if atom.name == 'MG  ' :
                MgIds.append((i_chain, i_res, i_atom))
            
'''    
######0         1         2         3         4         5         6         7         8
######012345678901234567890123456789012345678901234567890123456789012345678901234567890
   18 ATOM     18  P     G A   6     -17.691  53.155  69.007  1.00 78.49           P  
   19 ATOM     19  OP1   G A   6     -16.556  53.004  69.956  1.00 93.36           O  
   20 ATOM     20  OP2   G A   6     -18.664  52.043  68.845  1.00 95.58           O  
   21 ATOM     21  O5'   G A   6     -18.513  54.459  69.406  1.00 27.53           O 
   23 ATOM     23  C4'   G A   6     -18.227  56.897  69.444  1.00 20.08           C  
   24 ATOM     24  O4'   G A   6     -17.399  57.179  68.293  1.00 17.93           O  
   25 ATOM     25  C3'   G A   6     -19.643  57.145  68.949  1.00 17.69           C  
   26 ATOM     26  O3'   G A   6     -20.478  57.476  70.055  1.00 16.47           O  
   27 ATOM     27  C2'   G A   6     -19.468  58.342  68.013  1.00 17.00           C  
   28 ATOM     28  O2'   G A   6     -19.560  59.551  68.728  1.00 16.96           O  
   29 ATOM     29  C1'   G A   6     -18.022  58.183  67.527  1.00 16.40           C  
   30 ATOM     30  N9    G A   6     -17.823  57.884  66.116  1.00 48.19           N  
   31 ATOM     31  C8    G A   6     -16.924  58.497  65.277  1.00 47.26           C  
   32 ATOM     32  N7    G A   6     -16.970  58.042  64.054  1.00 46.20           N  
   33 ATOM     33  C5    G A   6     -17.959  57.071  64.087  1.00 46.43           C  
   34 ATOM     34  C6    G A   6     -18.462  56.251  63.055  1.00 46.33           C  
   35 ATOM     35  O6    G A   6     -18.136  56.229  61.865  1.00 49.61           O  
   36 ATOM     36  N1    G A   6     -19.451  55.396  63.525  1.00 45.10           N  
   37 ATOM     37  C2    G A   6     -19.900  55.342  64.826  1.00 46.79           C  
   38 ATOM     38  N2    G A   6     -20.853  54.434  65.091  1.00 43.52           N  
   39 ATOM     39  N3    G A   6     -19.448  56.118  65.796  1.00 49.01           N  
   40 ATOM     40  C4    G A   6     -18.486  56.952  65.359  1.00 47.62           C  
'''

def categorize_mp_RNA(atomname,res_name):
    if atomname[0] == 'H' :
        return 'hydrogen'
    
    elif atomname == ' O  ' and res_name == 'HOH' :
        return 'water_oxygen'
    
    elif atomname in ("MG  ","MN  ") :
        return 'metal'
    
    elif res_name in ("VAL", "LYS", "GLU", "LEU", "ALA", "GLY", "PHE", "HIS", "TYR", "ILE",
                      "ARG", "ASN", "THR", "MET", "ASP", "GLN", "TRP", "SER", "PRO", "HIS") :
        return 'protein'
    
    elif atomname in (" P  "," OP1"," OP2"," OP3") :
        return 'phosphate'
    
    elif atomname[3] == "'" :
        return 'sugar'
    
    else :
        return 'base'
                      
minset_base = []
minset_phos = []
minset_sugar= []
minset_water= []
minset_metal= []
minset_pro  = []

print (' the number of Mg : ',len(MgIds))
for icount, i in enumerate(MgIds) :
    print ('%i/%i' % (icount+1,len(MgIds)) )
    ''' 各Mgに対し、他の全残基に対する距離を計算'''
    #print (i)
    for (i_chain, chain) in enumerate(chains) :
        
        for (i_res, res) in enumerate(chain.residues) :
            '''各residueに対しては、各siteごとに最短の距離のみをとる'''
            '''複数のresidueにまたがって、一つのsiteを形成していることはない'''
            res_data = []
            '''とりあえず水素以外、residue内の全てのデータを入れて、あとでminを使う'''
            
            for (i_atom, atom) in enumerate(res.atoms) :
                if i == (i_chain, i_res, i_atom) : break
                
                category = categorize_mp_RNA(atom.name, atom.res_name)
                if category == 'hydrogen' : break
                
                distance = atom.xyz.distance(chains[i[0]].residues[i[1]].atoms[i[2]].xyz)
                #print (atom.name, atom.res_name, category, distance)
                res_data.append((category, distance))
                
            '''各site種類ごとに、このresidue内での最小値をminsetリストへ放りこむ'''
            if 'base' in [x[0] for x in res_data] :
                #print ('min(base)', min([x[1] for x in res_data if x[0] == 'base']))
                minset_base.append( min([x[1] for x in res_data if x[0] == 'base']) )
                
            if 'water_oxygen' in [x[0] for x in res_data]:
                #print ('min(water)', min([x[1] for x in res_data if x[0] == 'water_oxygen']))
                minset_water.append( min([x[1] for x in res_data if x[0] == 'water_oxygen']) )
                
            if 'sugar' in [x[0] for x in res_data] :
                #print ('min(sugar)', min([x[1] for x in res_data if x[0] == 'sugar']))
                minset_sugar.append( min([x[1] for x in res_data if x[0] == 'sugar']) )
                
            if 'phosphate' in [x[0] for x in res_data] :
                #print ('min(phos)', min([x[1] for x in res_data if x[0] == 'phosphate']))
                minset_phos.append( min([x[1] for x in res_data if x[0] == 'phosphate']) )
                
            if 'metal' in [x[0] for x in res_data] :
                #print ('min(metal)', min([x[1] for x in res_data if x[0] == 'metal']))
                minset_metal.append( min([x[1] for x in res_data if x[0] == 'metal']) )
                
            if 'protein' in [x[0] for x in res_data] :
                #print ('min(pro)', min([x[1] for x in res_data if x[0] == 'protein']))
                minset_pro.append( min([x[1] for x in res_data if x[0] == 'protein']) )
        
        
# gerenrate histograms
HIST_BIN = [HIST_MIN,]
x = HIST_MIN
i = 0
while x < HIST_MAX:
    i += 1 
    x = HIST_WID * i
    HIST_BIN.append(x)
#print (HIST_BIN)
hist_base = histogram(minset_base,  HIST_BIN)
hist_phos = histogram(minset_phos,  HIST_BIN)
hist_sugar= histogram(minset_sugar, HIST_BIN)
hist_water= histogram(minset_water, HIST_BIN)
hist_metal= histogram(minset_metal, HIST_BIN)
hist_pro  = histogram(minset_pro,   HIST_BIN)

# output to file_out
file_out.write('#0 Phosphate\n')
for i in xrange(len(hist_phos[0])) :
    file_out.write("%f %i\n" % (hist_phos[1][i], hist_phos[0][i]) )
file_out.write('\n\n#1 Base\n')
for i in xrange(len(hist_base[0])) :
    file_out.write("%f %i\n" % (hist_base[1][i], hist_base[0][i]) )
file_out.write('\n\n#2 Sugar\n')
for i in xrange(len(hist_sugar[0])) :
    file_out.write("%f %i\n" % (hist_sugar[1][i], hist_sugar[0][i]) )
file_out.write('\n\n#3 Oxygen of Water\n')
for i in xrange(len(hist_water[0])) :
    file_out.write("%f %i\n" % (hist_water[1][i], hist_water[0][i]) )
file_out.write('\n\n#4 Metal\n')
for i in xrange(len(hist_metal[0])) :
    file_out.write("%f %i\n" % (hist_metal[1][i], hist_metal[0][i]) )
file_out.write('\n\n#4 Protein\n')
for i in xrange(len(hist_pro[0])) :
    file_out.write("%f %i\n" % (hist_pro[1][i], hist_pro[0][i]) )
file_out.close()
                
