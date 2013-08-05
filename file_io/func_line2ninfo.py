#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@author: Naoto Hori
'''

from cafysis.elements.error import MyError
from cafysis.elements.ninfo import BondLength, BondAngle, Aicg13, Dihedral, AicgDih, Contact, BaseStack, BasePair

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

def line2aicg13(line) :
    it = iter(line.split())
    num = len(line.split())
    if num != 15 and num != 16:
        raise MyError("func_line2ninfo", "line2aicg13", "")
    if it.next() != 'aicg13' :
        raise MyError("func_line2ninfo", "line2aicg13", "This line is not aicg13.")
    info = Aicg13()
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
    info.wid = float(it.next())
    if num == 16:
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

def line2aicgdih(line) :
    it = iter(line.split())
    num = len(line.split())
    if num != 17 and num != 18:
        raise MyError("func_line2ninfo", "line2aicgdih", "")
    if it.next() != 'aicgdih' :
        raise MyError("func_line2ninfo", "line2aicgdih", "This line is not aicgdih.")
    info = AicgDih()
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
    info.wid = float(it.next())
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
