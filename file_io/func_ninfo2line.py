#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@author: Naoto Hori
'''

#0         1         2         3         4         5         6         7         8         9        10      
#012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567
#bond      1      1      1      1      2      1      2       4.8010       1.0000       1.0000      62.9100 SU
def bondlength2line(bl):
    s =  'bond'
    s += ' %6i' % (bl.id,)
    s += ' %6i %6i' % (bl.iunit1, bl.iunit2)
    s += ' %6i %6i' % (bl.imp1, bl.imp2)
    s += ' %6i %6i' % (bl.imp1un, bl.imp2un)
    s += ' %12.4f' % (bl.native,)
    s += ' %12.4f' % (bl.factor,)
    s += ' %12.4f' % (bl.correct_mgo,)
    s += ' %12.4f' % (bl.coef,)
    s += ' ' + bl.type
    s += '\n'
    return s 

#0         1         2         3         4         5         6         7         8         9        10      
#012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567
#fene      1      1      1      1      2      1      2       4.8010       2.0000       1.0000 SU
def fene2line(fene):
    s =  'fene'
    s += ' %6i' % (fene.id,)
    s += ' %6i %6i' % (fene.iunit1, fene.iunit2)
    s += ' %6i %6i' % (fene.imp1, fene.imp2)
    s += ' %6i %6i' % (fene.imp1un, fene.imp2un)
    s += ' %12.4f' % (fene.native,)
    s += ' %12.4f' % (fene.dist2,)
    s += ' %12.4f' % (fene.coef,)
    #s += ' ' + fene.type
    s += '\n'
    return s 

def bondangle2line(ba):
    s =  'angl'
    s += ' %6i' % (ba.id,)
    s += ' %6i %6i' % (ba.iunit1, ba.iunit2)
    s += ' %6i %6i %6i' % (ba.imp1, ba.imp2, ba.imp3)
    s += ' %6i %6i %6i' % (ba.imp1un, ba.imp2un, ba.imp3un)
    s += ' %12.4f' % (ba.native,)
    s += ' %12.4f' % (ba.factor,)
    s += ' %12.4f' % (ba.correct_mgo,)
    s += ' %12.4f' % (ba.coef,)
    s += ' ' + ba.type
    s += '\n'
    return s

#0         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        16
#01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
#aicg13      1      1      1      1      2      3      1      2      3       6.4310       1.0000       1.0000       1.1594       0.1500 ppp
def aicg132line(a13):
    s =  'aicg13'
    s += ' %6i' % (a13.id,)
    s += ' %6i %6i' % (a13.iunit1, a13.iunit2)
    s += ' %6i %6i %6i' % (a13.imp1, a13.imp2, a13.imp3)
    s += ' %6i %6i %6i' % (a13.imp1un, a13.imp2un, a13.imp3un)
    s += ' %12.4f' % (a13.native,)
    s += ' %12.4f' % (a13.factor,)
    s += ' %12.4f' % (a13.correct_mgo,)
    s += ' %12.4f' % (a13.coef,)
    s += ' %12.4f' % (a13.wid,)
    s += ' ' + a13.type
    s += '\n'
    return s 

def dihedral2line(dih):
    s =  'dihd'
    s += ' %6i' % (dih.id,)
    s += ' %6i %6i' % (dih.iunit1, dih.iunit2)
    s += ' %6i %6i %6i %6i' % (dih.imp1, dih.imp2, dih.imp3, dih.imp4)
    s += ' %6i %6i %6i %6i' % (dih.imp1un, dih.imp2un, dih.imp3un, dih.imp4un)
    s += ' %12.4f' % (dih.native,)
    s += ' %12.4f' % (dih.factor,)
    s += ' %12.4f' % (dih.correct_mgo,)
    s += ' %12.4f %12.4f' % (dih.coef, dih.coef_3)
    s += ' ' + dih.type
    s += '\n'
    return s 

def aicgdih2line(adih):
    s =  'aicgdih'
    s += ' %6i' % (adih.id,)
    s += ' %6i %6i' % (adih.iunit1, adih.iunit2)
    s += ' %6i %6i %6i %6i' % (adih.imp1, adih.imp2, adih.imp3, adih.imp4)
    s += ' %6i %6i %6i %6i' % (adih.imp1un, adih.imp2un, adih.imp3un, adih.imp4un)
    s += ' %12.4f' % (adih.native,)
    s += ' %12.4f' % (adih.factor,)
    s += ' %12.4f' % (adih.correct_mgo,)
    s += ' %12.4f %12.4f' % (adih.coef, adih.wid)
    s += ' ' + adih.type
    s += '\n'
    return s 

#**        icon iunit1-iunit2   imp1 - imp2 imp1un-imp2un      go_nat   factor_go  dummy     coef_go
#contact      5      1      1      2      5      2      5      6.1894      1.0000      1      0.9408 U-G
def contact2line(con):
    s =  'contact'
    s += ' %6i' % (con.id,)
    s += ' %6i %6i' % (con.iunit1, con.iunit2)
    s += ' %6i %6i' % (con.imp1, con.imp2)
    s += ' %6i %6i' % (con.imp1un, con.imp2un)
    s += ' %11.4f' % (con.native,)
    s += ' %11.4f' % (con.factor,)
    s += ' %6i' % (con.dummy,)
    s += ' %11.4f' % (con.coef,)
    s += ' ' + con.type
    s += '\n'
    return s 

#**     LJ iunit1-iunit2   imp1 - imp2 imp1un-imp2un    distance        coef
#LJ      5      1      1      2      5      2      5      6.1894      1.0000
def LJ2line(con):
    s =  'LJ'
    s += ' %6i' % (con.id,)
    s += ' %6i %6i' % (con.iunit1, con.iunit2)
    s += ' %6i %6i' % (con.imp1, con.imp2)
    s += ' %6i %6i' % (con.imp1un, con.imp2un)
    s += ' %11.4f' % (con.native,)
    s += ' %11.4f' % (con.coef,)
    s += '\n'
    return s 

def basepair2line(bp):
    s =  'basepair'
    s += ' %6i' % (bp.id,)
    s += ' %6i %6i' % (bp.iunit1, bp.iunit2)
    s += ' %6i %6i' % (bp.imp1, bp.imp2)
    s += ' %6i %6i' % (bp.imp1un, bp.imp2un)
    s += ' %11.4f' % (bp.native,)
    s += ' %11.4f' % (bp.factor,)
    s += ' %6i' % (bp.dummy,)
    s += ' %11.4f' % (bp.coef,)
    s += ' ' + bp.type
    s += ' %i' % (bp.nhb,)
    s += '\n'
    return s 

def basestack2line(bs):
    s =  'basestack'
    s += ' %6i' % (bs.id,)
    s += ' %6i %6i' % (bs.iunit1, bs.iunit2)
    s += ' %6i %6i' % (bs.imp1, bs.imp2)
    s += ' %6i %6i' % (bs.imp1un, bs.imp2un)
    s += ' %11.4f' % (bs.native,)
    s += ' %11.4f' % (bs.factor,)
    s += ' %6i' % (bs.dummy,)
    s += ' %11.4f' % (bs.coef,)
    s += ' ' + bs.type
    s += '\n'
    return s 
