#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2013/08/09
@author: Naoto Hori
'''

import sys
import Bio.PDB as PDB
    
TYPE2I = {'N':0, 'Ca':1, 'C':2, 'O':3, 'GCa':4, 'Cb':5,
          'KNz':6, 'KCd':7, 'DOd':8, 'RNh':9, 'NNd':10,
          'RNe':11, 'SOg':12, 'HNe':13, 'YCz':14, 'FCz':15,
          'LCd':16, 'CSg':17}

class TobiParam(object):
    def __init__(self,dirpath):
        self.r1 = {}
        self.r2 = {}
        self._read(dirpath+'/tobi_adp2_r1.ssv', 18, self.r1)
        self._check(self.r1)
        self._read(dirpath+'/tobi_adp2_r2.ssv', 18, self.r2)
        self._check(self.r2)
        
    def _read(self,filename, ntype, para):
        # top line
        f = open(filename, 'r')
        
        col = f.readline().split()[1:]
        for i in range(ntype):
            l = f.readline()
            lsp = l.split()
            if lsp[0] != col[i]:
                print('Error: lsp[0]!=col[0]')
                sys.exit(2)
            for m in range(ntype):
                para[(TYPE2I[col[i]],TYPE2I[col[m]])] = float(lsp[m+1])
        f.close()
        
    def _check(self,para):
        for mn in list(para.keys()):
            nm = (mn[1],mn[0])
            if para[mn] != para[nm]:
                print('Error: inconsistent')
                print('para[',mn,']=',para[mn])
                print('para[',nm,']=',para[nm])
    
def get_atom_type(res, atom):
    if atom == 'N':
        return TYPE2I['N']
    elif atom == 'CA':
        if res == 'GLY':
            return TYPE2I['GCa']
        else:
            return TYPE2I['Ca']
    elif atom == 'C':
        return TYPE2I['C']
    elif atom == 'O':
        return TYPE2I['O']
    elif ((atom == 'CB' and res != 'SER') or
          (atom in ('CG','CD') and res == 'PRO')):
        return TYPE2I['Cb']
    elif res == 'LYS' and atom in ('CE','NZ'):
        return TYPE2I['KNz']
    elif res == 'LYS' and atom == 'CD':
        return TYPE2I['KCd']
    elif ((res == 'ASP') or
          (res == 'GLU' and atom in ('CD','OE1','OE2'))):
        return TYPE2I['DOd']
    elif res == 'ARG' and atom in ('CZ','NH1','NH2'):
        return TYPE2I['RNh']
    elif ((res == 'ASN')or
          (res == 'GLN' and atom in ('CD','OE1','NE2'))):
        return TYPE2I['NNd']
    elif res == 'ARG' and atom in ('CD','NE'):
        return TYPE2I['RNe']
    elif ((res == 'SER' and atom in ('CB','OG')) or
          (res == 'THR' and atom == 'OG1') or
          (res == 'TYR' and atom == 'OH')):
        return TYPE2I['SOg']
    elif ((res == 'HIS' and atom in ('CG','ND1','CD2','CE1','NE2')) or
          (res == 'TRP' and atom == 'NE1')):
        return TYPE2I['HNe']
    elif res == 'TYR' and atom in ('CE1', 'CE2','CZ'):
        return TYPE2I['YCz']
    elif ((atom == 'CG' and res in ('ARG','GLN','GLU','LEU','LYS','MET','PHE','TRP','TYR')) or 
          (atom == 'CG1' and res == 'ILE') or
          (atom == 'SD' and res == 'MET') or
          (res == 'PHE' and atom in ('CD1','CD2','CE1','CE2','CZ')) or
          (res == 'THR' and atom == 'CG2') or
          (res == 'TRP' and atom in ('CD1','CD2','CE2','CE3','CZ2','CZ3','CH2')) or
          (res == 'TYR' and atom in ('CD1','CD2'))):
        return TYPE2I['FCz']
    elif ((res == 'ILE' and atom in ('CG2','CD1')) or 
          (res == 'LEU' and atom in ('CD1','CD2')) or
          (res == 'MET' and atom == 'CE') or
          (res == 'VAL' and atom in ('CG1','CG2'))):
        return TYPE2I['LCd']
    elif res == 'CYS' and atom == 'SG':
        return TYPE2I['CSg']
    else:
        print('Error: get_atom_type in tobi.py')
        sys.exit(2)
        

if __name__ == "__main__":
    
    
    if len(sys.argv) != 2:
        print('Usage: SCRIPT [input PDB]')
        sys.exit(2)
        
    import os
    #params = TobiParam('~/python/cafysis/param/')
    params = TobiParam(os.path.dirname(__file__) + '/para/')
        
    p = PDB.PDBParser(PERMISSIVE=1)
    ## PDB読み込みのデバッグではPERMISSIVE=0でも試す
    
    #structs = []
    #for i,filename in enumerate(sys.argv[1:]):
    #    structs.append(p.get_structure(i,filename))
    
    #struct = p.get_structure('complex',sys.argv[1])
    chains = p.get_structure('complex',sys.argv[1])[0].get_list()
        
    def accept_residue(res):
        'アミノ酸のみTrue'
        #ResidueのID(タプル)の第一要素(hetero)が空白以外
        if res.get_id()[0] != ' ': 
            return False
        if res.get_resname() in ('ARG','GLY','ALA','ASP','ASN',
                                 'GLU','GLN','PHE','LYS','ILE',
                                 'LEU','TRP','CYS','THR','HIS',
                                 'MET','PRO','VAL','TYR','SER',
                                 'HIE'):
            return True
        else:
            return False
        
    import re
    re_H = re.compile("^[\dH]")  # 数字またはHではじまる場合は、水素と判定
    def accept_atom(atom):
        '水素以外True'
        aname = atom.get_name()
        if re_H.match(aname):
            return False
        else:
            return True
        
    
    energy = {}
    num_chain = len(chains)
            
    atoms_chain = []
    types_chain = []
    for c in chains:
        atoms = []
        types = []
        for a in c.get_atoms():
            r = a.get_parent()
            rname = r.get_resname()
            aname = a.get_name()
            if not (accept_residue(r) and accept_atom(a)):
                #print '############', resname1, aname1
                continue
            atoms.append(a)
            types.append(get_atom_type(rname, aname))
        atoms_chain.append(atoms)
        types_chain.append(types)
    
    print('start main')
    
    for i in range(num_chain):
        for j in range(i+1, num_chain):
            energy[(i,j)] = 0.0
            
            #ii=0
            for a1,t1 in zip(atoms_chain[i],types_chain[i]):
                #ii+=1
                #jj=0
                for a2,t2 in zip(atoms_chain[j],types_chain[j]):
                    #jj+=1
                    dist = a2 - a1
                    if dist > 6.0:
                        #ene = 0.0
                        pass
                    elif dist <= 4.0:
                        #ene = params.r1[(t1,t2)]
                        #energy[(i,j)] += ene
                        energy[(i,j)] += params.r1[(t1,t2)]
                    elif dist <= 6.0:
                        #ene = params.r2[(t1,t2)]
                        #energy[(i,j)] += ene
                        energy[(i,j)] += params.r2[(t1,t2)]
                    #print '%i %i %i %i %10.5f' % (ii, jj, t1, t2, ene)
                    
    for i in range(num_chain):
        for j in range(i+1,num_chain):
            print('%i %i %12.6f' % (i,j,energy[(i,j)]))
                    
            
