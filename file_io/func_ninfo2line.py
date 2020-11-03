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
    #s += ' %12.4f' % (fene.dist2,)
    s += ' %12.4f' % (fene.factor,)
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

def basestackDTdist2line(bs):
    s =  'bs-dist'
    s += ' %6i' % (bs.id,)
    s += ' %6i %6i' % (bs.iunit1, bs.iunit2)
    s += ' %6i %6i' % (bs.imp1, bs.imp2)
    s += ' %6i %6i' % (bs.imp1un, bs.imp2un)
    s += ' %8.4f' % (bs.h,)
    s += ' %8.4f' % (bs.s,)
    s += ' %6.2f' % (bs.Tm,)
    s += ' %11.4f' % (bs.native,)
    s += ' %11.4f' % (bs.coef,)
    s += ' ' + bs.type
    s += '\n'
    return s 

def basestackDTdih2line(bs):
    s =  'bs-dihd'
    s += ' %6i' % (bs.id,)
    s += ' %6i' % (bs.dih1_id,)
    s += ' %6i %6i' % (bs.dih1_iunit1, bs.dih1_iunit2)
    s += ' %6i %6i %6i %6i' % (bs.dih1_imp1, bs.dih1_imp2, bs.dih1_imp3, bs.dih1_imp4)
    s += ' %6i %6i %6i %6i' % (bs.dih1_imp1un, bs.dih1_imp2un, bs.dih1_imp3un, bs.dih1_imp4un)
    s += ' %11.4f' % (bs.dih1_native,)
    s += ' %11.4f' % (bs.dih1_coef,)
    s += ' ' + bs.dih1_type
    s += '\n'
    s += 'bs-dihd'
    s += ' %6i' % (bs.id,)
    s += ' %6i' % (bs.dih2_id,)
    s += ' %6i %6i' % (bs.dih2_iunit1, bs.dih2_iunit2)
    s += ' %6i %6i %6i %6i' % (bs.dih2_imp1, bs.dih2_imp2, bs.dih2_imp3, bs.dih2_imp4)
    s += ' %6i %6i %6i %6i' % (bs.dih2_imp1un, bs.dih2_imp2un, bs.dih2_imp3un, bs.dih2_imp4un)
    s += ' %11.4f' % (bs.dih2_native,)
    s += ' %11.4f' % (bs.dih2_coef,)
    s += ' ' + bs.dih2_type
    s += '\n'
    return s 

def hbondDTdist2line(hb):
    s =  'hb-dist'
    s += ' %6i' % (hb.id,)
    s += ' %6i %6i' % (hb.iunit1, hb.iunit2)
    s += ' %6i %6i' % (hb.imp1, hb.imp2)
    s += ' %6i %6i' % (hb.imp1un, hb.imp2un)
    s += ' %11.4f' % (hb.factor,)
    s += ' %11.4f' % (hb.native,)
    s += ' %11.4f' % (hb.coef,)
    if hb.sectert:
        s += ' %s' % hb.sectert
    if hb.nHB:
        s += ' %i' % hb.nHB
    if hb.atoms1 and hb.atoms1:
        for a1, a2 in zip(hb.atoms1, hb.atoms2):
            s += ' %s %s' % (a1,a2)
    s += '\n'
    return s 

def hbondDTangl2line(hb):
    s =  'hb-angl'
    s += ' %6i' % (hb.id,)
    s += ' %6i' % (hb.ang1_id,)
    s += ' %6i %6i' % (hb.ang1_iunit1, hb.ang1_iunit2)
    s += ' %6i %6i %6i' % (hb.ang1_imp1, hb.ang1_imp2, hb.ang1_imp3)
    s += ' %6i %6i %6i' % (hb.ang1_imp1un, hb.ang1_imp2un, hb.ang1_imp3un)
    s += ' %11.4f' % (hb.ang1_native,)
    s += ' %11.4f' % (hb.ang1_coef,)
    s += '\n'
    s += 'hb-angl'
    s += ' %6i' % (hb.id,)
    s += ' %6i' % (hb.ang2_id,)
    s += ' %6i %6i' % (hb.ang2_iunit1, hb.ang2_iunit2)
    s += ' %6i %6i %6i' % (hb.ang2_imp1, hb.ang2_imp2, hb.ang2_imp3)
    s += ' %6i %6i %6i' % (hb.ang2_imp1un, hb.ang2_imp2un, hb.ang2_imp3un)
    s += ' %11.4f' % (hb.ang2_native,)
    s += ' %11.4f' % (hb.ang2_coef,)
    s += '\n'
    return s 

def hbondDTdih2line(hb):
    s =  'hb-dihd'
    s += ' %6i' % (hb.id,)
    s += ' %6i' % (hb.dih0_id,)
    s += ' %6i %6i' % (hb.dih0_iunit1, hb.dih0_iunit2)
    s += ' %6i %6i %6i %6i' % (hb.dih0_imp1, hb.dih0_imp2, hb.dih0_imp3, hb.dih0_imp4)
    s += ' %6i %6i %6i %6i' % (hb.dih0_imp1un, hb.dih0_imp2un, hb.dih0_imp3un, hb.dih0_imp4un)
    s += ' %11.4f' % (hb.dih0_native,)
    s += ' %11.4f' % (hb.dih0_coef,)
    s += '\n'
    s += 'hb-dihd'
    s += ' %6i' % (hb.id,)
    s += ' %6i' % (hb.dih1_id,)
    s += ' %6i %6i' % (hb.dih1_iunit1, hb.dih1_iunit2)
    s += ' %6i %6i %6i %6i' % (hb.dih1_imp1, hb.dih1_imp2, hb.dih1_imp3, hb.dih1_imp4)
    s += ' %6i %6i %6i %6i' % (hb.dih1_imp1un, hb.dih1_imp2un, hb.dih1_imp3un, hb.dih1_imp4un)
    s += ' %11.4f' % (hb.dih1_native,)
    s += ' %11.4f' % (hb.dih1_coef,)
    s += '\n'
    s += 'hb-dihd'
    s += ' %6i' % (hb.id,)
    s += ' %6i' % (hb.dih2_id,)
    s += ' %6i %6i' % (hb.dih2_iunit1, hb.dih2_iunit2)
    s += ' %6i %6i %6i %6i' % (hb.dih2_imp1, hb.dih2_imp2, hb.dih2_imp3, hb.dih2_imp4)
    s += ' %6i %6i %6i %6i' % (hb.dih2_imp1un, hb.dih2_imp2un, hb.dih2_imp3un, hb.dih2_imp4un)
    s += ' %11.4f' % (hb.dih2_native,)
    s += ' %11.4f' % (hb.dih2_coef,)
    s += '\n'
    return s 

def tertiarystackDTdist2line(tst):
    s =  'tbs-dist'
    s += ' %6i' % (tst.id,)
    s += ' %6i %6i' % (tst.iunit1, tst.iunit2)
    s += ' %6i %6i' % (tst.imp1, tst.imp2)
    s += ' %6i %6i' % (tst.imp1un, tst.imp2un)
    s += ' %11.4f' % (tst.factor,)
    s += ' %11.4f' % (tst.native,)
    s += ' %11.4f' % (tst.coef,)
    s += ' %+2i %+2i' % (tst.excess1, tst.excess2)
    s += '\n'
    return s 

def tertiarystackDTangl2line(tst):
    s =  'tbs-angl'
    s += ' %6i' % (tst.id,)
    s += ' %6i' % (tst.ang1_id,)
    s += ' %6i %6i' % (tst.ang1_iunit1, tst.ang1_iunit2)
    s += ' %6i %6i %6i' % (tst.ang1_imp1, tst.ang1_imp2, tst.ang1_imp3)
    s += ' %6i %6i %6i' % (tst.ang1_imp1un, tst.ang1_imp2un, tst.ang1_imp3un)
    s += ' %11.4f' % (tst.ang1_native,)
    s += ' %11.4f' % (tst.ang1_coef,)
    s += '\n'
    s += 'tbs-angl'
    s += ' %6i' % (tst.id,)
    s += ' %6i' % (tst.ang2_id,)
    s += ' %6i %6i' % (tst.ang2_iunit1, tst.ang2_iunit2)
    s += ' %6i %6i %6i' % (tst.ang2_imp1, tst.ang2_imp2, tst.ang2_imp3)
    s += ' %6i %6i %6i' % (tst.ang2_imp1un, tst.ang2_imp2un, tst.ang2_imp3un)
    s += ' %11.4f' % (tst.ang2_native,)
    s += ' %11.4f' % (tst.ang2_coef,)
    s += '\n'
    return s 

def tertiarystackDTdih2line(tst):
    s =  'tbs-dihd'
    s += ' %6i' % (tst.id,)
    s += ' %6i' % (tst.dih0_id,)
    s += ' %6i %6i' % (tst.dih0_iunit1, tst.dih0_iunit2)
    s += ' %6i %6i %6i %6i' % (tst.dih0_imp1, tst.dih0_imp2, tst.dih0_imp3, tst.dih0_imp4)
    s += ' %6i %6i %6i %6i' % (tst.dih0_imp1un, tst.dih0_imp2un, tst.dih0_imp3un, tst.dih0_imp4un)
    s += ' %11.4f' % (tst.dih0_native,)
    s += ' %11.4f' % (tst.dih0_coef,)
    s += '\n'
    s += 'tbs-dihd'
    s += ' %6i' % (tst.id,)
    s += ' %6i' % (tst.dih1_id,)
    s += ' %6i %6i' % (tst.dih1_iunit1, tst.dih1_iunit2)
    s += ' %6i %6i %6i %6i' % (tst.dih1_imp1, tst.dih1_imp2, tst.dih1_imp3, tst.dih1_imp4)
    s += ' %6i %6i %6i %6i' % (tst.dih1_imp1un, tst.dih1_imp2un, tst.dih1_imp3un, tst.dih1_imp4un)
    s += ' %11.4f' % (tst.dih1_native,)
    s += ' %11.4f' % (tst.dih1_coef,)
    s += '\n'
    s += 'tbs-dihd'
    s += ' %6i' % (tst.id,)
    s += ' %6i' % (tst.dih2_id,)
    s += ' %6i %6i' % (tst.dih2_iunit1, tst.dih2_iunit2)
    s += ' %6i %6i %6i %6i' % (tst.dih2_imp1, tst.dih2_imp2, tst.dih2_imp3, tst.dih2_imp4)
    s += ' %6i %6i %6i %6i' % (tst.dih2_imp1un, tst.dih2_imp2un, tst.dih2_imp3un, tst.dih2_imp4un)
    s += ' %11.4f' % (tst.dih2_native,)
    s += ' %11.4f' % (tst.dih2_coef,)
    s += '\n'
    return s 
