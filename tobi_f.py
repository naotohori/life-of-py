#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2013/08/09
@author: Naoto Hori
'''

import sys
import Bio.PDB as PDB
import numpy as np
import py_calc_tobi
    
TYPE2I = {'N':0, 'Ca':1, 'C':2, 'O':3, 'GCa':4, 'Cb':5,
          'KNz':6, 'KCd':7, 'DOd':8, 'RNh':9, 'NNd':10,
          'RNe':11, 'SOg':12, 'HNe':13, 'YCz':14, 'FCz':15,
          'LCd':16, 'CSg':17}

class TobiParam(object):
    def __init__(self,dirpath):
        self.r1 = {}
        self.r2 = {}
        self.NUM_PARA = 18
        self._read(dirpath+'/tobi_adp2_r1.ssv', self.NUM_PARA, self.r1)
        self._check(self.r1)
        self._read(dirpath+'/tobi_adp2_r2.ssv', self.NUM_PARA, self.r2)
        self._check(self.r2)
        self.param1 = np.zeros((self.NUM_PARA, self.NUM_PARA))
        self.param2 = np.zeros((self.NUM_PARA, self.NUM_PARA))
        for i in xrange(self.NUM_PARA):
            for j in xrange(self.NUM_PARA):
                self.param1[i,j] = self.r1[(i,j)]
                self.param2[i,j] = self.r2[(i,j)]
        
    def _read(self,filename, ntype, para):
        # top line
        f = open(filename, 'r')
        
        col = f.readline().split()[1:]
        for i in xrange(ntype):
            l = f.readline()
            lsp = l.split()
            if lsp[0] != col[i]:
                print 'Error: lsp[0]!=col[0]'
                sys.exit(2)
            for m in xrange(ntype):
                para[(TYPE2I[col[i]],TYPE2I[col[m]])] = float(lsp[m+1])
        f.close()
        
    def _check(self,para):
        for mn in para.keys():
            nm = (mn[1],mn[0])
            if para[mn] != para[nm]:
                print 'Error: inconsistent'
                print 'para[',mn,']=',para[mn]
                print 'para[',nm,']=',para[nm]
    
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
        print 'Error: get_atom_type in tobi.py'
        sys.exit(2)
        
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
        
def accept_atom(atom):
    import re
    re_H = re.compile("^[\dH]")  # 数字またはHではじまる場合は、水素と判定
    '水素以外True'
    aname = atom.get_name()
    if re_H.match(aname):
        return False
    else:
        return True
        
    
def calc_tobi_for_pdb(pdb_filepath):
    p = PDB.PDBParser(PERMISSIVE=1)
    ## PDB読み込みのデバッグではPERMISSIVE=0でも試す
    
    #structs = []
    #for i,filename in enumerate(sys.argv[1:]):
    #    structs.append(p.get_structure(i,filename))
    
    #struct = p.get_structure('complex',sys.argv[1])
    chains = p.get_structure('complex',pdb_filepath)[0].get_list()
    
    num_chain = len(chains)
            
    '''accept or not'''
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
    
    '''list to ndarray'''
    num_atom = [len(x) for x in atoms_chain]
    xyz = np.zeros((3,max(num_atom),num_chain))
    atom2type = np.zeros( (max(num_atom),num_chain), dtype=np.int)
    
    for ichain in xrange(num_chain):
        for iatom in xrange(num_atom[ichain]):
            c = atoms_chain[ichain][iatom].get_coord()
            xyz[0,iatom,ichain] = c[0]
            xyz[1,iatom,ichain] = c[1]
            xyz[2,iatom,ichain] = c[2]
            atom2type[iatom,ichain] = types_chain[ichain][iatom]
        
    #print 'start main'
    '''call the Fortran module'''
    ene = py_calc_tobi.calc_tobi( params.param1, params.param2, 4.0, 6.0,
                                  num_atom, atom2type, xyz)
    
    return ene

if __name__ == "__main__":
    
    if not len(sys.argv) in (2,3):
        print 'Usage: SCRIPT [input PDB]'
        print 'Usage: SCRIPT [input PDB dir] [output file]'
        sys.exit(2)
        
    import os
    
    #params = TobiParam('~/python/cafysis/param/')
    params = TobiParam(os.path.dirname(__file__) + '/para/')
        
    if len(sys.argv) == 2:
        ene = calc_tobi_for_pdb(sys.argv[1])
        
        for i in xrange(len(ene)):
            for j in xrange(i+1,len(ene)):
                print i+1,j+1,ene[i,j]    
                
    elif len(sys.argv) == 3:
        import glob
        outfile = open(sys.argv[2], 'w')
        
        pdbfiles = glob.glob(sys.argv[1]+'*.pdb')
        pdbfiles.sort()
        for pdbfilepath in pdbfiles:
            
            ene = calc_tobi_for_pdb(pdbfilepath)
            
            for i in xrange(len(ene)):
                for j in xrange(i+1,len(ene)):
                    outfile.write('%s %i %i %12.6f\n' 
                      % ( os.path.basename(pdbfilepath)[:-4], i+1, j+1, ene[i,j]))
            
        outfile.close()
            