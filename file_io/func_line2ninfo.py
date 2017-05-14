#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@author: Naoto Hori
'''

from cafysis.elements.error import MyError
from cafysis.elements.ninfo import BondLength, Fene, BondAngle, Dihedral, Contact, LJ, BaseStack, BasePair, BaseStackDT, TertiaryStackDT, HBondDT

def line2bondlength(line) :
    it = iter(line.split())
    num = len(line.split())
    if num != 12 and num != 13:
        raise MyError("func_line2ninfo", "line2bondlength", "")
    if it.next() != 'bond' :
        raise MyError("func_line2ninfo", "line2bondlength", "This line is not bondlength.")
    info = BondLength()
    info.id = int(it.next())
    info.iunit1 = int(it.next())
    info.iunit2 = int(it.next())
    info.imp1 = int(it.next())
    info.imp2 = int(it.next())
    info.imp1un = int(it.next())
    info.imp2un = int(it.next())
    info.native = float(it.next())
    info.factor = float(it.next())
    info.correct_mgo = float(it.next())
    info.coef = float(it.next())
    if num == 13:
        info.type = it.next()
    return info

def line2bondangle(line) :
    it = iter(line.split())
    num = len(line.split())
    if num != 14 and num != 15:
        raise MyError("func_line2ninfo", "line2bondangle", "")
    if it.next() != 'angl' :
        raise MyError("func_line2ninfo", "line2bondangle", "This line is not bondangle.")
    info = BondAngle()
    info.id = int(it.next())
    info.iunit1 = int(it.next())
    info.iunit2 = int(it.next())
    info.imp1 = int(it.next())
    info.imp2 = int(it.next())
    info.imp3 = int(it.next())
    info.imp1un = int(it.next())
    info.imp2un = int(it.next())
    info.imp3un = int(it.next())
    info.native = float(it.next())
    info.factor = float(it.next())
    info.correct_mgo = float(it.next())
    info.coef = float(it.next())
    if num == 15:
        info.type = it.next()
    return info

 
def line2dihedral(line) :
    it = iter(line.split())
    num = len(line.split())
    if num != 17 and num != 18:
        raise MyError("func_line2ninfo", "line2dihedral", "")
    if it.next() != 'dihd' :
        raise MyError("func_line2ninfo", "line2dihedral", "This line is not dihedral.")
    info = Dihedral()
    info.id = int(it.next())
    info.iunit1 = int(it.next())
    info.iunit2 = int(it.next())
    info.imp1 = int(it.next())
    info.imp2 = int(it.next())
    info.imp3 = int(it.next())
    info.imp4 = int(it.next())
    info.imp1un = int(it.next())
    info.imp2un = int(it.next())
    info.imp3un = int(it.next())
    info.imp4un = int(it.next())
    info.native = float(it.next())
    info.factor = float(it.next())
    info.correct_mgo = float(it.next())
    info.coef = float(it.next())
    info.coef_3 = float(it.next())
    if num == 18:
        info.type = it.next()
    return info


def line2contact(line) :
    it = iter(line.split())
    num = len(line.split())
    if num != 12 and num != 13:
        raise MyError("func_line2ninfo", "line2contact", "")
    if it.next() != 'contact' :
        raise MyError("func_line2ninfo", "line2contact", "This line is not contact.")
    info = Contact()
    info.id = int(it.next())
    info.iunit1 = int(it.next())
    info.iunit2 = int(it.next())
    info.imp1 = int(it.next())
    info.imp2 = int(it.next())
    info.imp1un = int(it.next())
    info.imp2un = int(it.next())
    info.native = float(it.next())
    info.factor = float(it.next())
    info.dummy = int(it.next())
    info.coef = float(it.next())
    if num == 13:
        info.type = it.next()
    return info

def line2basepair(line) :
    it = iter(line.split())
    num = len(line.split())
    if num != 14 :
        raise MyError("func_line2ninfo", "line2basepair", "")
    if it.next() != 'basepair' :
        raise MyError("func_line2ninfo", "line2basepair", "This line is not basepair")
    info = BasePair()
    info.id = int(it.next())
    info.iunit1 = int(it.next())
    info.iunit2 = int(it.next())
    info.imp1 = int(it.next())
    info.imp2 = int(it.next())
    info.imp1un = int(it.next())
    info.imp2un = int(it.next())
    info.native = float(it.next())
    info.factor = float(it.next())
    info.dummy = int(it.next())
    info.coef = float(it.next())
    info.type = it.next()
    info.nhb = int(it.next())
    return info

def line2basestack(line) :
    it = iter(line.split())
    num = len(line.split())
    if num != 13 :
        raise MyError("func_line2ninfo", "line2basestack", "")
    if it.next() != 'basestack' :
        raise MyError("func_line2ninfo", "line2basestack", "This line is not basestack")
    info = BaseStack()
    info.id = int(it.next())
    info.iunit1 = int(it.next())
    info.iunit2 = int(it.next())
    info.imp1 = int(it.next())
    info.imp2 = int(it.next())
    info.imp1un = int(it.next())
    info.imp2un = int(it.next())
    info.native = float(it.next())
    info.factor = float(it.next())
    info.dummy = int(it.next())
    info.coef = float(it.next())
    info.type = it.next()
    return info


def line2basestackDT(line) :
    it = iter(line.split())
    num = len(line.split())
    if num != 12 :
        raise MyError("func_line2ninfo", "line2basestackDT", "")
    if it.next() != 'bs-dist':
        raise MyError("func_line2ninfo", "line2basestackDT", "This line is not bs-dist")
    info = BaseStackDT()
    info.id = int(it.next())
    info.iunit1 = int(it.next())
    info.iunit2 = int(it.next())
    info.imp1 = int(it.next())
    info.imp2 = int(it.next())
    info.imp1un = int(it.next())
    info.imp2un = int(it.next())
    info.factor = float(it.next())
    info.native = float(it.next())
    info.coef = float(it.next())
    info.type = it.next()
    info.dih1_id = None
    info.dih2_id = None
    return info


def line2tertiarystackDT(line) :
    it = iter(line.split())
    num = len(line.split())
    if num != 13 :
        raise MyError("func_line2ninfo", "line2tertiarystackDT", "")
    if it.next() != 'tbs-dist':
        raise MyError("func_line2ninfo", "line2tertiarystackDT", "This line is not tbs-dist")
    info = TertiaryStackDT()
    info.id = int(it.next())
    info.iunit1 = int(it.next())
    info.iunit2 = int(it.next())
    info.imp1 = int(it.next())
    info.imp2 = int(it.next())
    info.imp1un = int(it.next())
    info.imp2un = int(it.next())
    info.factor = float(it.next())
    info.native = float(it.next())
    info.coef = float(it.next())
    info.excess1 = int(it.next())
    info.excess2 = int(it.next())

    info.ang1_id = None
    info.ang2_id = None
    info.dih0_id = None
    info.dih1_id = None
    info.dih2_id = None
    return info


def line2hbondDT(line) :
    it = iter(line.split())
    num = len(line.split())
    if num < 14 :
        raise MyError("func_line2ninfo", "line2hbondDT", "")
    if it.next() != 'hb-dist':
        raise MyError("func_line2ninfo", "line2hbondDT", "This line is not hb-dist")
    info = HBondDT()
    info.id = int(it.next())
    info.iunit1 = int(it.next())
    info.iunit2 = int(it.next())
    info.imp1 = int(it.next())
    info.imp2 = int(it.next())
    info.imp1un = int(it.next())
    info.imp2un = int(it.next())
    info.factor = float(it.next())
    info.native = float(it.next())
    info.coef = float(it.next())
    info.sectert = it.next()
    info.nHB = int(it.next())
    info.atoms1 = []
    info.atoms2 = []
    try:
        while (True):
            a = it.next()
            info.atoms1.append(a)
            a = it.next()
            info.atoms2.append(a)
    except:
        pass

    info.ang1_id = None
    info.ang2_id = None
    info.dih0_id = None
    info.dih1_id = None
    info.dih2_id = None
    return info
