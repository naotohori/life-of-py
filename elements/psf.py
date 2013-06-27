#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2011/12/08
@author: Naoto Hori
'''

class Atom(object):
    def __init__(self):
        self.atom_id = None
        self.seg_name = None
        self.res_id = None
        self.res_name = None
        self.atom_name = None
        self.atom_type = None
        self.charge = 0.0
        self.mass = 0.0
        self.unused = 0
    
    def show(self):
        print (self.atom_id, self.seg_name, 
               self.res_id, self.res_name,
               self.atom_name, self.atom_type,
               self.charge, self.mass, self.unused)

class Psf(object):
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.atoms = []
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.impropers = []
        self.crossterms = []
        self.n_atom = lambda : len(self.atoms)
        self.n_bond = lambda : len(self.bonds)
        self.n_angle = lambda : len(self.angles)
        self.n_dihedral = lambda : len(self.dihedrals)
        self.n_improper = lambda : len(self.impropers)
        self.n_crossterm = lambda : len(self.crossterms)
        
    def show(self):
        for atom in self.atoms :
            atom.show()
        print self.bonds
        print self.angles
        print self.dihedrals
        print self.impropers
        print self.crossterms
        