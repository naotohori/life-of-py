#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@author: Naoto Hori
'''

from cafysis.elements.error import MyError
from cafysis.file_io.func_line2ninfo import *
from cafysis.file_io.func_ninfo2line import *

class NinfoFile(object):
    def __init__(self, filename):
        self._filename = filename
        self._status = 'Closed'
        
    def open_to_read(self):
        if self._status != 'Closed' :
            raise MyError('NinfoFile', 'open_for_read', 'file is not closed')
        self._file = open(self._filename, 'r')
        
    def open_to_write(self):
        if self._status != 'Closed' :
            raise MyError('NinfoFile', 'open_for_read', 'file is not closed')
        self._file = open(self._filename, 'w')

    def close(self):
        self._file.close()
        
    def read_all(self, ni):
        for line in self._file :
            if line[0:4] == 'bond' :
                ni.bondlengths.append(line2bondlength(line))
            elif line[0:4] == 'fene' :
                ni.fenes.append(line2fene(line))
            elif line[0:4] == 'angl' :
                ni.bondangles.append(line2bondangle(line))
            elif line[0:4] == 'dihd' :
                ni.dihedrals.append(line2dihedral(line))
            elif line[0:7] == 'contact' :
                ni.contacts.append(line2contact(line))
            elif line[0:2] == 'LJ' :
                ni.contacts.append(line2LJ(line))
            elif line[0:9] == 'basestack' :
                ni.basestacks.append(line2basestack(line))
            elif line[0:8] == 'basepair' :
                ni.basepairs.append(line2basepair(line))

    def write_all(self, ni):
        ni.update_info()
        num_unit = ni.max_unit
        
        if len(ni.bondlengths) > 0:
            self._file.write('<<<< native bond length\n')
            self._file.write('**      ibd iunit1-iunit2   imp1 - imp2 imp1un-imp2un      bd_nat    factor_bd  correct_mgo      coef_bd\n')
            for bl in ni.bondlengths :
                self._file.write(bondlength2line(bl))
            self._file.write('>>>>\n')
            self._file.write('\n')

        if len(ni.fenes) > 0:
            self._file.write('<<<< FENE\n')
            self._file.write('**      ibd iunit1-iunit2   imp1 - imp2 imp1un-imp2un        dist        dist2        coef\n')
            for fene in ni.fenes :
                self._file.write(fene2line(fene))
            self._file.write('>>>>\n')
            self._file.write('\n')
        
        if len(ni.bondangles) > 0:
            self._file.write('<<<< native bond angles\n')
            self._file.write('**      iba iunit1-iunit2   imp1 - imp2 - imp3 imp1un-imp2un-imp3un      ba_nat    factor_ba  correct_mgo      coef_ba\n')
            for ba in ni.bondangles :
                self._file.write(bondangle2line(ba))
            self._file.write('>>>>\n')
            self._file.write('\n')
        
        if len(ni.dihedrals) > 0:
            self._file.write('<<<< native dihedral angles\n')
            self._file.write('**     idih iunit1-iunit2   imp1 - imp2 - imp3 - imp4 imp1un-imp2un-imp3un-imp4un      dih_nat   factor_dih  correct_mgo   coef_dih_1   coef_dih_3\n')
            for dih in ni.dihedrals :
                self._file.write(dihedral2line(dih))
            self._file.write('>>>>\n')
            self._file.write('\n')
            
        if len(ni.contacts) > 0:
            self._file.write('<<<< native contact\n')
            self._file.write('** total_contact = %i\n'%(len(ni.contacts),))
            for i in xrange(1, num_unit+1) :
                for j in xrange(i, num_unit+1) :
                    subset = ni.get_contacts_by_unit(i,j)
                    if len(subset) == 0:
                        continue
                    self._file.write('\n')
                    self._file.write('** contact between unit %6i and %6i\n'%(i,j))
                    self._file.write('** total_contact_unit = %i\n'%(len(subset),))
                    self._file.write('**        icon iunit1-iunit2   imp1 - imp2 imp1un-imp2un      go_nat   factor_go  dummy     coef_go\n')
                    for con in subset:
                        self._file.write(contact2line(con))
            self._file.write('>>>>\n')
            self._file.write('\n')

        if len(ni.LJs) > 0:
            self._file.write('<<<< LJ\n')
            self._file.write('** total_LJ = %i\n'%(len(ni.LJs),))
            for i in xrange(1, num_unit+1) :
                for j in xrange(i, num_unit+1) :
                    subset = ni.get_LJs_by_unit(i,j)
                    if len(subset) == 0:
                        continue
                    self._file.write('\n')
                    self._file.write('** LJ between unit %6i and %6i\n'%(i,j))
                    self._file.write('** total_LJ_unit = %i\n'%(len(subset),))
                    self._file.write('**     iLJ  iunit1-iunit2   imp1 - imp2 imp1un-imp2un  distance        coef\n')
                    for con in subset:
                        self._file.write(LJ2line(con))
            self._file.write('>>>>\n')
            self._file.write('\n')
        
        if len(ni.basepairs) > 0:
            self._file.write('<<<< native basepair\n')
            self._file.write('** total_contact = %i\n'%(len(ni.basepairs),))
            for i in xrange(1, num_unit+1) :
                for j in xrange(i, num_unit+1) :
                    subset = ni.get_basepairs_by_unit(i,j)
                    if len(subset) == 0:
                        continue
                    self._file.write('\n')
                    self._file.write('** contact between unit %6i and %6i\n'%(i,j))
                    self._file.write('** total_contact_unit = %i\n'%(len(subset),))
                    self._file.write('**        icon iunit1-iunit2   imp1 - imp2 imp1un-imp2un      go_nat   factor_go  dummy     coef_go\n')
                    for bp in subset:
                        self._file.write(basepair2line(bp))
            self._file.write('>>>>\n')
            self._file.write('\n')
        
        if len(ni.basestacks) > 0:
            self._file.write('<<<< native basestack\n')
            self._file.write('** total_contact = %i\n'%(len(ni.basestacks),))
            self._file.write('**        icon iunit1-iunit2   imp1 - imp2 imp1un-imp2un      go_nat   factor_go  dummy     coef_go\n')
            for bs in ni.basestacks :
                self._file.write(basestack2line(bs))
            self._file.write('>>>>\n')
        
    def write_unit(self, ni, un1, un2):
        self._file.write('<<<< native bond length\n')
        self._file.write('**      ibd iunit1-iunit2   imp1 - imp2 imp1un-imp2un      bd_nat    factor_bd  correct_mgo      coef_bd\n')
        for bl in ni.bondlengths :
            if bl.iunit1 == un1 and bl.iunit2 == un2:
                self._file.write(bondlength2line(bl))
        self._file.write('>>>>\n')
        self._file.write('\n')

        self._file.write('<<<< FENE\n')
        self._file.write('**      ibd iunit1-iunit2   imp1 - imp2 imp1un-imp2un      bd_nat    factor_bd  correct_mgo      coef_bd\n')
        for fene in ni.fenes :
            if fene.iunit1 == un1 and fene.iunit2 == un2:
                self._file.write(fene2line(fene))
        self._file.write('>>>>\n')
        self._file.write('\n')
        
        self._file.write('<<<< native bond angles\n')
        self._file.write('**      iba iunit1-iunit2   imp1 - imp2 - imp3 imp1un-imp2un-imp3un      ba_nat    factor_ba  correct_mgo      coef_ba\n')
        for ba in ni.bondangles :
            if ba.iunit1 == un1 and ba.iunit2 == un2:
                self._file.write(bondangle2line(ba))
        self._file.write('>>>>\n')
        self._file.write('\n')
        
        self._file.write('<<<< native dihedral angles\n')
        self._file.write('**     idih iunit1-iunit2   imp1 - imp2 - imp3 - imp4 imp1un-imp2un-imp3un-imp4un      dih_nat   factor_dih  correct_mgo   coef_dih_1   coef_dih_3\n')
        for dih in ni.dihedrals :
            if dih.iunit1 == un1 and dih.iunit2 == un2:
                self._file.write(dihedral2line(dih))
        self._file.write('>>>>\n')
        self._file.write('\n')
            
        ni.update_info()
        
        subset = ni.get_contacts_by_unit(un1,un2)
        self._file.write('<<<< native contact\n')
        self._file.write('** contact between unit %6i and %6i\n'%(un1,un2))
        self._file.write('** total_contact_unit = %i\n'%(len(subset),))
        self._file.write('**        icon iunit1-iunit2   imp1 - imp2 imp1un-imp2un      go_nat   factor_go  dummy     coef_go\n')
        for con in subset:
            self._file.write(contact2line(con))
        self._file.write('>>>>\n')
        self._file.write('\n')


        subset = ni.get_LJs_by_unit(un1,un2)
        self._file.write('<<<< LJ\n')
        self._file.write('** LJ between unit %6i and %6i\n'%(un1,un2))
        self._file.write('** total_LJ_unit = %i\n'%(len(subset),))
        self._file.write('**     iLJ  iunit1-iunit2   imp1 - imp2 imp1un-imp2un    distance        coef\n')
        for con in subset:
            self._file.write(LJ2line(con))
        self._file.write('>>>>\n')
        self._file.write('\n')
        
        subset = ni.get_basepairs_by_unit(un1,un2)
        if len(subset) != 0:
            self._file.write('<<<< native basepair\n')
            self._file.write('\n')
            self._file.write('** contact between unit %6i and %6i\n'%(un1,un2))
            self._file.write('** total_contact_unit = %i\n'%(len(subset),))
            self._file.write('**        icon iunit1-iunit2   imp1 - imp2 imp1un-imp2un      go_nat   factor_go  dummy     coef_go\n')
            for bp in subset:
                self._file.write(basepair2line(bp))
            self._file.write('>>>>\n')
            self._file.write('\n')
        
        self._file.write('<<<< native basestack\n')
        self._file.write('** total_contact = %i\n'%(len(ni.basestacks),))
        self._file.write('**        icon iunit1-iunit2   imp1 - imp2 imp1un-imp2un      go_nat   factor_go  dummy     coef_go\n')
        for bs in ni.basestacks :
            if bs.iunit1==un1 and bs.iunit2 == un2:
                self._file.write(basestack2line(bs))
        self._file.write('>>>>\n')
