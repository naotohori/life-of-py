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
    if next(it) != 'bond' :
        raise MyError("func_line2ninfo", "line2bondlength", "This line is not bondlength.")
    info = BondLength()
    info.id = int(next(it))
    info.iunit1 = int(next(it))
    info.iunit2 = int(next(it))
    info.imp1 = int(next(it))
    info.imp2 = int(next(it))
    info.imp1un = int(next(it))
    info.imp2un = int(next(it))
    info.native = float(next(it))
    info.factor = float(next(it))
    info.correct_mgo = float(next(it))
    info.coef = float(next(it))
    if num == 13:
        info.type = next(it)
    return info

def line2bondangle(line) :
    it = iter(line.split())
    num = len(line.split())
    if num != 14 and num != 15:
        raise MyError("func_line2ninfo", "line2bondangle", "")
    if next(it) != 'angl' :
        raise MyError("func_line2ninfo", "line2bondangle", "This line is not bondangle.")
    info = BondAngle()
    info.id = int(next(it))
    info.iunit1 = int(next(it))
    info.iunit2 = int(next(it))
    info.imp1 = int(next(it))
    info.imp2 = int(next(it))
    info.imp3 = int(next(it))
    info.imp1un = int(next(it))
    info.imp2un = int(next(it))
    info.imp3un = int(next(it))
    info.native = float(next(it))
    info.factor = float(next(it))
    info.correct_mgo = float(next(it))
    info.coef = float(next(it))
    if num == 15:
        info.type = next(it)
    return info

 
def line2dihedral(line) :
    it = iter(line.split())
    num = len(line.split())
    if num != 17 and num != 18:
        raise MyError("func_line2ninfo", "line2dihedral", "")
    if next(it) != 'dihd' :
        raise MyError("func_line2ninfo", "line2dihedral", "This line is not dihedral.")
    info = Dihedral()
    info.id = int(next(it))
    info.iunit1 = int(next(it))
    info.iunit2 = int(next(it))
    info.imp1 = int(next(it))
    info.imp2 = int(next(it))
    info.imp3 = int(next(it))
    info.imp4 = int(next(it))
    info.imp1un = int(next(it))
    info.imp2un = int(next(it))
    info.imp3un = int(next(it))
    info.imp4un = int(next(it))
    info.native = float(next(it))
    info.factor = float(next(it))
    info.correct_mgo = float(next(it))
    info.coef = float(next(it))
    info.coef_3 = float(next(it))
    if num == 18:
        info.type = next(it)
    return info


def line2contact(line) :
    it = iter(line.split())
    num = len(line.split())
    if num != 12 and num != 13:
        raise MyError("func_line2ninfo", "line2contact", "")
    if next(it) != 'contact' :
        raise MyError("func_line2ninfo", "line2contact", "This line is not contact.")
    info = Contact()
    info.id = int(next(it))
    info.iunit1 = int(next(it))
    info.iunit2 = int(next(it))
    info.imp1 = int(next(it))
    info.imp2 = int(next(it))
    info.imp1un = int(next(it))
    info.imp2un = int(next(it))
    info.native = float(next(it))
    info.factor = float(next(it))
    info.dummy = int(next(it))
    info.coef = float(next(it))
    if num == 13:
        info.type = next(it)
    return info

def line2basepair(line) :
    it = iter(line.split())
    num = len(line.split())
    if num != 14 :
        raise MyError("func_line2ninfo", "line2basepair", "")
    if next(it) != 'basepair' :
        raise MyError("func_line2ninfo", "line2basepair", "This line is not basepair")
    info = BasePair()
    info.id = int(next(it))
    info.iunit1 = int(next(it))
    info.iunit2 = int(next(it))
    info.imp1 = int(next(it))
    info.imp2 = int(next(it))
    info.imp1un = int(next(it))
    info.imp2un = int(next(it))
    info.native = float(next(it))
    info.factor = float(next(it))
    info.dummy = int(next(it))
    info.coef = float(next(it))
    info.type = next(it)
    info.nhb = int(next(it))
    return info

def line2basestack(line) :
    it = iter(line.split())
    num = len(line.split())
    if num != 13 :
        raise MyError("func_line2ninfo", "line2basestack", "")
    if next(it) != 'basestack' :
        raise MyError("func_line2ninfo", "line2basestack", "This line is not basestack")
    info = BaseStack()
    info.id = int(next(it))
    info.iunit1 = int(next(it))
    info.iunit2 = int(next(it))
    info.imp1 = int(next(it))
    info.imp2 = int(next(it))
    info.imp1un = int(next(it))
    info.imp2un = int(next(it))
    info.native = float(next(it))
    info.factor = float(next(it))
    info.dummy = int(next(it))
    info.coef = float(next(it))
    info.type = next(it)
    return info


def line2basestackDT(line, fmt=2) :

    if fmt == 1:
        it = iter(line.split())
        num = len(line.split())
        if num != 12 :
            raise MyError("func_line2ninfo", "line2basestackDT", "")
        if next(it) != 'bs-dist':
            raise MyError("func_line2ninfo", "line2basestackDT", "This line is not bs-dist")
        info = BaseStackDT()
        info.id = int(next(it))
        info.iunit1 = int(next(it))
        info.iunit2 = int(next(it))
        info.imp1 = int(next(it))
        info.imp2 = int(next(it))
        info.imp1un = int(next(it))
        info.imp2un = int(next(it))
        info.factor = float(next(it))
        info.native = float(next(it))
        info.coef = float(next(it))
        info.type = next(it)
        info.dih1_id = None
        info.dih2_id = None
        return info

    elif fmt == 2:
        it = iter(line.split())
        num = len(line.split())
        if num != 14 :
            raise MyError("func_line2ninfo", "line2basestackDT", "num!=14")
        if next(it) != 'bs-dist':
            raise MyError("func_line2ninfo", "line2basestackDT", "This line is not bs-dist")
        info = BaseStackDT()
        info.id = int(next(it))
        info.iunit1 = int(next(it))
        info.iunit2 = int(next(it))
        info.imp1 = int(next(it))
        info.imp2 = int(next(it))
        info.imp1un = int(next(it))
        info.imp2un = int(next(it))
        info.h = float(next(it))
        info.s = float(next(it))
        info.Tm = float(next(it))
        info.native = float(next(it))
        info.coef = float(next(it))
        info.type = next(it)
        info.dih1_id = None
        info.dih2_id = None
        return info

    else:
        raise MyError("func_line2ninfo", "line2basestackDT", "unknown fmt")


def line2tertiarystackDT(line) :
    it = iter(line.split())
    num = len(line.split())
    if num != 13 :
        raise MyError("func_line2ninfo", "line2tertiarystackDT", "")
    if next(it) != 'tbs-dist':
        raise MyError("func_line2ninfo", "line2tertiarystackDT", "This line is not tbs-dist")
    info = TertiaryStackDT()
    info.id = int(next(it))
    info.iunit1 = int(next(it))
    info.iunit2 = int(next(it))
    info.imp1 = int(next(it))
    info.imp2 = int(next(it))
    info.imp1un = int(next(it))
    info.imp2un = int(next(it))
    info.factor = float(next(it))
    info.native = float(next(it))
    info.coef = float(next(it))
    info.excess1 = int(next(it))
    info.excess2 = int(next(it))

    info.ang1_id = None
    info.ang2_id = None
    info.dih0_id = None
    info.dih1_id = None
    info.dih2_id = None
    return info


def line2hbondDT(line) :
    it = iter(line.split())
    num = len(line.split())
    if num < 11 :
        raise MyError("func_line2ninfo", "line2hbondDT", "")
    if next(it) != 'hb-dist':
        raise MyError("func_line2ninfo", "line2hbondDT", "This line is not hb-dist")
    info = HBondDT()
    info.id = int(next(it))
    info.iunit1 = int(next(it))
    info.iunit2 = int(next(it))
    info.imp1 = int(next(it))
    info.imp2 = int(next(it))
    info.imp1un = int(next(it))
    info.imp2un = int(next(it))
    info.factor = float(next(it))
    info.native = float(next(it))
    info.coef = float(next(it))
    if num > 11:
        info.sectert = next(it)
        info.nHB = int(next(it))
        info.atoms1 = []
        info.atoms2 = []
        try:
            while (True):
                a = next(it)
                info.atoms1.append(a)
                a = next(it)
                info.atoms2.append(a)
        except:
            pass
    else:
        info.sectert = None
        info.nHB = None

    info.ang1_id = None
    info.ang2_id = None
    info.dih0_id = None
    info.dih1_id = None
    info.dih2_id = None
    return info
