#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2011/12/08
@author: Naoto Hori
'''
from lop.elements.error import MyError
from lop.elements.psf import Psf, Atom

def line2atom(line):
    line.rstrip()
    atom = Atom()
    ## Old
    #atom.atom_id = int(line[0:8])
    #atom.seg_name = line[8:11]
    #atom.res_id = int(line[11:15])
    #atom.res_name = line[15:23]
    #atom.atom_name = line[23:28]
    #atom.atom_type = line[28:31]
    #atom.charge = float(line[31:50])
    #atom.mass = float(line[50:64])
    #atom.unused = int(line[64:76])
    #atom.atom_id = int(line[0:8])
    #atom.seg_name = line[8:10]
    #atom.res_id = int(line[10:15])
    #atom.res_name = line[15:23]
    #atom.atom_name = line[23:29]
    #atom.atom_type = line[29:33]
    #atom.charge = float(line[33:51])
    #atom.mass = float(line[51:64])
    #atom.unused = int(line[64:76])
    ## New
    atom.atom_id = int(line[0:8])
    atom.seg_name = line[8:13]
    atom.res_id = int(line[13:18])
    atom.res_name = line[18:23]
    atom.atom_name = line[22:27]
    atom.atom_type = line[27:31]
    atom.charge = float(line[31:44])
    atom.mass = float(line[44:58])
    atom.unused = int(line[58:70])
    return atom

def atom2line(atom):
    ## Old
    return ('%8i%3s%4i%8s%5s%3s%19.6f%14.4f%12i' %
    ## New
    #return ('%8i  %3s %4i  %3s %4s %4s  %9.6f      %8.4f%12i' %
             (atom.atom_id, atom.seg_name,
              atom.res_id, atom.res_name,
              atom.atom_name, atom.atom_type,
              atom.charge, atom.mass,
              atom.unused) )
    
class PsfFile(object):
    '''
    classdocs
    '''

    def __init__(self, filename):
        '''
        Constructor
        '''
        self._filename = filename
        self._status = 'Closed'
        
    def open_to_read(self):
        if self._status != 'Closed' :
            raise MyError('PsfFile', 'open_for_read', 'file is not closed')
        self._file = open(self._filename, 'r')
        
    def open_to_write(self):
        if self._status != 'Closed' :
            raise MyError('PsfFile', 'open_for_read', 'file is not closed')
        self._file = open(self._filename, 'w')
        
    def close(self):
        self._file.close()
        
    def read_all(self):
        psf = Psf()
        
        header = self._file.readline()
        while header :
            if header.find('!NATOM') != -1 :
                n = int(header.split()[0])
                for i in range(n) :
                    line = self._file.readline()
                    psf.atoms.append(line2atom(line))
            if header.find('!NBOND') != -1 :
                n = int(header.split()[0])
                i = 0
                while i < n :
                    a = [int(x) for x in self._file.readline().split()]
                    for j in range(0, len(a), 2) :
                        psf.bonds.append((a[j],a[j+1]))
                    i += len(a) / 2
            if header.find('!NTHETA') != -1 :
                n = int(header.split()[0])
                i = 0
                while i < n :
                    a = [int(x) for x in self._file.readline().split()]
                    for j in range(0, len(a), 3) :
                        psf.angles.append((a[j],a[j+1],a[j+2]))
                    i += len(a) / 3
            if header.find('!NPHI') != -1 :
                n = int(header.split()[0])
                i = 0
                while i < n :
                    a = [int(x) for x in self._file.readline().split()]
                    for j in range(0, len(a), 4) :
                        psf.dihedrals.append((a[j],a[j+1],a[j+2],a[j+3]))
                    i += len(a) / 4
            if header.find('!NIMPHI') != -1 :
                n = int(header.split()[0])
                i = 0
                while i < n :
                    a = [int(x) for x in self._file.readline().split()]
                    for j in range(0, len(a), 4) :
                        psf.impropers.append((a[j],a[j+1],a[j+2],a[j+3]))
                    i += len(a) / 4
            if header.find('!NCRTERM') != -1 :
                n = int(header.split()[0])
                i = 0
                while i < n :
                    a = [int(x) for x in self._file.readline().split()]
                    for j in range(0, len(a), 4) :
                        psf.crossterms.append((a[j],a[j+1],a[j+2],a[j+3]))
                    i += len(a) / 4
            header = self._file.readline()
        
        return psf
    
    def write_all(self, psf):
        if psf.n_atom() > 0:
            self._file.write('%8i !NATOM\n' % psf.n_atom())
            for atom in psf.atoms :
                self._file.write(atom2line(atom) + '\n')
        if psf.n_bond() > 0 :
            self._file.write('%8i !NBOND\n' % psf.n_bond())
            for i in range(0, psf.n_bond()) :
                self._file.write('%8i%8i' % psf.bonds[i])
                if (i+1)%4==0 or i==psf.n_bond()-1 :
                    self._file.write('\n')
        if psf.n_angle() > 0 :
            self._file.write('%8i !NTHETA\n' % psf.n_angle())
            for i in range(0, psf.n_angle()) :
                self._file.write('%8i%8i%8i' % psf.angles[i])
                if (i+1)%3==0 or i==psf.n_angle()-1 :
                    self._file.write('\n')
        if psf.n_dihedral() > 0 :
            self._file.write('%8i !NPHI\n' % psf.n_dihedral())
            for i in range(0, psf.n_dihedral()) :
                self._file.write('%8i%8i%8i' % psf.dihedrals[i])
                if (i+1)%2==0 or i==psf.n_dihedral()-1 :
                    self._file.write('\n')
        if psf.n_improper() > 0 :
            self._file.write('%8i !NIMPHI\n' % psf.n_improper())
            for i in range(0, psf.n_improper()) :
                self._file.write('%8i%8i%8i' % psf.impropers[i])
                if (i+1)%2==0 or i==psf.n_improper()-1 :
                    self._file.write('\n')
        if psf.n_crossterm() > 0 :
            self._file.write('%8i !NCRTERM\n' % psf.n_crossterm())
            for i in range(0, psf.n_crossterm()) :
                self._file.write('%8i%8i%8i' % psf.crossterms[i])
                if (i+1)%2==0 or i==psf.n_crossterm()-1 :
                    self._file.write('\n')
    
if __name__ == '__main__':
    import sys
    
    if len(sys.argv) != 3:
        print('This is debug modeg mode of PsfFile class')
        print('Usage: SCRIPT [input psf file] [output psf file]')
        sys.exit(2)
    
    in_file = PsfFile(sys.argv[1])    
    in_file.open_to_read()
    psf = in_file.read_all()
    psf.show()
    out_file = PsfFile(sys.argv[2])
    out_file.open_to_write()
    out_file.write_all(psf)
