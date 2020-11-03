#!/usr/bin/env python
'''
@author: Naoto Hori
'''
import copy

class Ninfo(object):
    def __init__(self,id=None,iunit1=None,iunit2=None, imp1=None,imp2=None,imp1un=None,imp2un=None,
                 native=None,factor=None,correct_mgo=None,coef=None,type_str=None):
        self.id = id
        self.iunit1 = iunit1
        self.iunit2 = iunit2
        self.imp1 = imp1
        self.imp2 = imp2
        self.imp1un = imp1un
        self.imp2un = imp2un
        self.native = native
        self.factor = factor
        self.correct_mgo = correct_mgo
        self.coef = coef
        self.type = type_str

    def show(self):
        print('id',self.id)
        print('iunit1',self.iunit1)
        print('iunit2',self.iunit2)
        print('imp1',self.imp1)
        print('imp2',self.imp2)
        print('imp1un',self.imp1un)
        print('imp2un',self.imp2un)
        print('native',self.native)
        print('factor',self.factor)
        print('correct_mgo',self.correct_mgo)
        print('coef',self.coef)
        print('type',self.type)
        
#**      ibd iunit1-iunit2   imp1 - imp2 imp1un-imp2un      bd_nat    factor_bd  correct_mgo      coef_bd
#bond      1      1      1      1      2      1      2       3.8531       1.0000       1.0000     110.4000 pp
class BondLength(Ninfo):
    _counter = 0
    def __init__(self,id=None,iunit1=None,iunit2=None, imp1=None,imp2=None,imp1un=None,imp2un=None,
                 native=None,factor=None,correct_mgo=None,coef=None,type_str=None):
        if id is None:
            BondLength._counter += 1
            id = BondLength._counter
        Ninfo.__init__(self,id=id,iunit1=iunit1,iunit2=iunit2, imp1=imp1,imp2=imp2,imp1un=imp1un,imp2un=imp2un,
                       native=native,factor=factor,correct_mgo=correct_mgo,coef=coef,type_str=type_str)
            
        
#**      ibd iunit1-iunit2   imp1 - imp2 imp1un-imp2un      bd_nat         dist2         coef
#fene      1      1      1      1      2      1      2       3.8531       2.0000       1.0000
class Fene(Ninfo):
    _counter = 0
    def __init__(self,id=None,iunit1=None,iunit2=None, imp1=None,imp2=None,imp1un=None,imp2un=None,
                 native=None,factor=None,correct_mgo=None,coef=None,type_str=None):
        if id is None:
            Fene._counter += 1
            id = Fene._counter
        Ninfo.__init__(self,id=id,iunit1=iunit1,iunit2=iunit2, imp1=imp1,imp2=imp2,imp1un=imp1un,imp2un=imp2un,
                       native=native,factor=factor,correct_mgo=correct_mgo,coef=coef,type_str=type_str)
        self.dist2 = None

    
class BondAngle(Ninfo):
#**      iba iunit1-iunit2   imp1 - imp2 - imp3 imp1un-imp2un-imp3un      ba_nat    factor_ba  correct_mgo      coef_ba
#angl      1      1      1      1      2      3      1      2      3     114.0257       0.0000      1.0000       0.0000 ppp
    _counter = 0
    def __init__(self,id=None,iunit1=None,iunit2=None, imp1=None,imp2=None,imp1un=None,imp2un=None,imp3=None,imp3un=None,
                 native=None,factor=None,correct_mgo=None,coef=None,type_str=None):
        if id is None:
            BondAngle._counter += 1
            id = BondAngle._counter
        Ninfo.__init__(self,id=id,iunit1=iunit1,iunit2=iunit2, imp1=imp1,imp2=imp2,imp1un=imp1un,imp2un=imp2un,
                       native=native,factor=factor,correct_mgo=correct_mgo,coef=coef,type_str=type_str)
        self.imp3 = imp3
        self.imp3un = imp3un
    
    def show(self):
        Ninfo.show(self)
        print('imp3',self.imp3)
        print('imp3un',self.imp3un)

class Dihedral(Ninfo):
#**     idih iunit1-iunit2   imp1 - imp2 - imp3 - imp4 imp1un-imp2un-imp3un-imp4un      dih_nat   factor_dih  correct_mgo   coef_dih_1   coef_dih_3
#dihd      1      1      1      1      2      3      4      1      2      3      4    -132.8457       0.0000       1.0000       0.0000       0.0000 pppp
    _counter = 0
    def __init__(self,id=None,iunit1=None,iunit2=None, imp1=None,imp2=None,imp1un=None,imp2un=None,
                 imp3=None,imp3un=None,imp4=None,imp4un=None,coef_3=None,
                 native=None,factor=None,correct_mgo=None,coef=None,type_str=None):
        if id is None:
            Dihedral._counter += 1
            id = Dihedral._counter
        Ninfo.__init__(self,id=id,iunit1=iunit1,iunit2=iunit2, imp1=imp1,imp2=imp2,imp1un=imp1un,imp2un=imp2un,
                       native=native,factor=factor,correct_mgo=correct_mgo,coef=coef,type_str=type_str)
        self.imp3 = imp3
        self.imp4 = imp4
        self.imp3un = imp3un
        self.imp4un = imp4un
        self.coef_3 = coef_3
    
    def show(self):
        Ninfo.show(self)
        print('imp3',self.imp3)
        print('imp4',self.imp4)
        print('imp3un',self.imp3un)
        print('imp4un',self.imp4un)
        print('coef_3',self.coef_3)

class Contact(Ninfo):
    _counter = 0
    def __init__(self,id=None,iunit1=None,iunit2=None, imp1=None,imp2=None,imp1un=None,imp2un=None,dummy=None,
                 native=None,factor=None,correct_mgo=None,coef=None,type_str=None):
        if id is None:
            Contact._counter += 1
            id = Contact._counter
        Ninfo.__init__(self,id=id,iunit1=iunit1,iunit2=iunit2, imp1=imp1,imp2=imp2,imp1un=imp1un,imp2un=imp2un,
                       native=native,factor=factor,correct_mgo=correct_mgo,coef=coef,type_str=type_str)
        self.dummy = dummy

class LJ(Ninfo):
    _counter = 0 
    def __init__(self,id=None,iunit1=None,iunit2=None, imp1=None,imp2=None,imp1un=None,imp2un=None,
                 native=None,factor=None,correct_mgo=None,coef=None,type_str=None):
        if id is None:
            LJ._counter += 1
            id = LJ._counter
        Ninfo.__init__(self,id=id,iunit1=iunit1,iunit2=iunit2, imp1=imp1,imp2=imp2,imp1un=imp1un,imp2un=imp2un,
                       native=native,factor=factor,correct_mgo=correct_mgo,coef=coef,type_str=type_str)
    
class BasePair(Ninfo):
    _counter = 0
    def __init__(self,id=None,iunit1=None,iunit2=None, imp1=None,imp2=None,imp1un=None,imp2un=None,dummy=None,nhb=None,
                 native=None,factor=None,correct_mgo=None,coef=None,type_str=None):
        if id is None:
            BasePair._counter += 1
            id = BasePair._counter
        Ninfo.__init__(self,id=id,iunit1=iunit1,iunit2=iunit2, imp1=imp1,imp2=imp2,imp1un=imp1un,imp2un=imp2un,
                       native=native,factor=factor,correct_mgo=correct_mgo,coef=coef,type_str=type_str)
        self.dummy = dummy
        self.nhb = nhb
    
class BaseStack(Ninfo):
    _counter = 0
    def __init__(self,id=None,iunit1=None,iunit2=None, imp1=None,imp2=None,imp1un=None,imp2un=None,dummy=None,
                 native=None,factor=None,correct_mgo=None,coef=None,type_str=None):
        if id is None:
            BaseStack._counter += 1
            id = BaseStack._counter
        Ninfo.__init__(self,id=id,iunit1=iunit1,iunit2=iunit2, imp1=imp1,imp2=imp2,imp1un=imp1un,imp2un=imp2un,
                       native=native,factor=factor,correct_mgo=correct_mgo,coef=coef,type_str=type_str)
        self.dummy = dummy

class BaseStackDT(Ninfo):
    _counter = 0
    _counter_dih = 0
    def __init__(self,id=None,id_dih1=None,id_dih2=None,
                 iunit1=None,iunit2=None, imp1=None,imp2=None,imp1un=None,imp2un=None,
                 h=None, s=None, Tm=None,
                 native=None,factor=None,correct_mgo=None,coef=None,type_str=None,
                 dih1_imp1=None, dih1_imp2=None, dih1_imp3=None, dih1_imp4=None,dih1_iunit1=None,dih1_iunit2=None,
                 dih1_imp1un=None, dih1_imp2un=None, dih1_imp3un=None, dih1_imp4un=None,
                 dih1_native=None,dih1_coef=None,dih1_type_str=None,
                 dih2_imp1=None, dih2_imp2=None, dih2_imp3=None, dih2_imp4=None,dih2_iunit1=None,dih2_iunit2=None,
                 dih2_imp1un=None, dih2_imp2un=None, dih2_imp3un=None, dih2_imp4un=None,
                 dih2_native=None,dih2_coef=None,dih2_type_str=None):

        if id is None:
            BaseStackDT._counter += 1
            id = BaseStackDT._counter
        if id_dih1 is None:
            BaseStackDT._counter_dih += 1
            id_dih1 = BaseStackDT._counter_dih
        if id_dih2 is None:
            BaseStackDT._counter_dih += 1
            id_dih2 = BaseStackDT._counter_dih

        Ninfo.__init__(self,id=id,iunit1=iunit1,iunit2=iunit2, imp1=imp1,imp2=imp2,imp1un=imp1un,imp2un=imp2un,
                       native=native,factor=factor,correct_mgo=correct_mgo,coef=coef,type_str=type_str)
        self.h = h
        self.s = s
        self.Tm = Tm
        self.dih1_id = id_dih1
        self.dih1_iunit1 = dih1_iunit1
        self.dih1_iunit2 = dih1_iunit2
        self.dih1_imp1 = dih1_imp1
        self.dih1_imp2 = dih1_imp2
        self.dih1_imp3 = dih1_imp3
        self.dih1_imp4 = dih1_imp4
        self.dih1_imp1un = dih1_imp1un
        self.dih1_imp2un = dih1_imp2un
        self.dih1_imp3un = dih1_imp3un
        self.dih1_imp4un = dih1_imp4un
        self.dih1_native = dih1_native
        self.dih1_coef = dih1_coef
        self.dih1_type = dih1_type_str
        self.dih2_id = id_dih2
        self.dih2_iunit1 = dih2_iunit1
        self.dih2_iunit2 = dih2_iunit2
        self.dih2_imp1 = dih2_imp1
        self.dih2_imp2 = dih2_imp2
        self.dih2_imp3 = dih2_imp3
        self.dih2_imp4 = dih2_imp4
        self.dih2_imp1un = dih2_imp1un
        self.dih2_imp2un = dih2_imp2un
        self.dih2_imp3un = dih2_imp3un
        self.dih2_imp4un = dih2_imp4un
        self.dih2_native = dih2_native
        self.dih2_coef = dih2_coef
        self.dih2_type = dih2_type_str
    
class HBondDT(Ninfo):
    _counter = 0
    _counter_ang = 0
    _counter_dih = 0
    def __init__(self,id=None,id_ang1=None,id_ang2=None,id_dih0=None,id_dih1=None,id_dih2=None,
                 iunit1=None,iunit2=None, imp1=None,imp2=None,imp1un=None,imp2un=None,
                 native=None,factor=None,correct_mgo=None,coef=None,type_str=None,
                 ang1_iunit1=None,ang1_iunit2=None,
                 ang1_imp1=None,ang1_imp2=None,ang1_imp1un=None,ang1_imp2un=None,ang1_imp3=None,ang1_imp3un=None,
                 ang1_native=None,ang1_coef=None,ang1_type_str=None,
                 ang2_iunit1=None,ang2_iunit2=None,
                 ang2_imp1=None,ang2_imp2=None,ang2_imp1un=None,ang2_imp2un=None,ang2_imp3=None,ang2_imp3un=None,
                 ang2_native=None,ang2_coef=None,ang2_type_str=None,
                 dih0_imp1=None, dih0_imp2=None, dih0_imp3=None, dih0_imp4=None,dih0_iunit1=None,dih0_iunit2=None,
                 dih0_imp1un=None, dih0_imp2un=None, dih0_imp3un=None, dih0_imp4un=None,
                 dih0_native=None,dih0_coef=None,dih0_type_str=None,
                 dih1_imp1=None, dih1_imp2=None, dih1_imp3=None, dih1_imp4=None,dih1_iunit1=None,dih1_iunit2=None,
                 dih1_imp1un=None, dih1_imp2un=None, dih1_imp3un=None, dih1_imp4un=None,
                 dih1_native=None,dih1_coef=None,dih1_type_str=None,
                 dih2_imp1=None, dih2_imp2=None, dih2_imp3=None, dih2_imp4=None,dih2_iunit1=None,dih2_iunit2=None,
                 dih2_imp1un=None, dih2_imp2un=None, dih2_imp3un=None, dih2_imp4un=None,
                 dih2_native=None,dih2_coef=None,dih2_type_str=None, sectert=None, nHB=None, atoms1=None, atoms2=None):
        if id is None:
            HBondDT._counter += 1
            id = HBondDT._counter
        if id_ang1 is None:
            HBondDT._counter_ang += 1
            id_ang1 = HBondDT._counter_ang
        if id_ang2 is None:
            HBondDT._counter_ang += 1
            id_ang2 = HBondDT._counter_ang
        if id_dih0 is None:
            HBondDT._counter_dih += 1
            id_dih0 = HBondDT._counter_dih
        if id_dih1 is None:
            HBondDT._counter_dih += 1
            id_dih1 = HBondDT._counter_dih
        if id_dih2 is None:
            HBondDT._counter_dih += 1
            id_dih2 = HBondDT._counter_dih
        Ninfo.__init__(self,id=id,iunit1=iunit1,iunit2=iunit2, imp1=imp1,imp2=imp2,imp1un=imp1un,imp2un=imp2un,
                       native=native,factor=factor,correct_mgo=correct_mgo,coef=coef,type_str=type_str)
        self.ang1_id = id_ang1
        self.ang1_iunit1 = ang1_iunit1
        self.ang1_iunit2 = ang1_iunit2
        self.ang1_imp1 = ang1_imp1
        self.ang1_imp2 = ang1_imp2
        self.ang1_imp3 = ang1_imp3
        self.ang1_imp1un = ang1_imp1un
        self.ang1_imp2un = ang1_imp2un
        self.ang1_imp3un = ang1_imp3un
        self.ang1_native = ang1_native
        self.ang1_coef = ang1_coef
        self.ang1_type = ang1_type_str
        self.ang2_id = id_ang2
        self.ang2_iunit1 = ang2_iunit1
        self.ang2_iunit2 = ang2_iunit2
        self.ang2_imp1 = ang2_imp1
        self.ang2_imp2 = ang2_imp2
        self.ang2_imp3 = ang2_imp3
        self.ang2_imp1un = ang2_imp1un
        self.ang2_imp2un = ang2_imp2un
        self.ang2_imp3un = ang2_imp3un
        self.ang2_native = ang2_native
        self.ang2_coef = ang2_coef
        self.ang2_type = ang2_type_str
        self.dih0_id = id_dih0
        self.dih0_iunit1 = dih0_iunit1
        self.dih0_iunit2 = dih0_iunit2
        self.dih0_imp1 = dih0_imp1
        self.dih0_imp2 = dih0_imp2
        self.dih0_imp3 = dih0_imp3
        self.dih0_imp4 = dih0_imp4
        self.dih0_imp1un = dih0_imp1un
        self.dih0_imp2un = dih0_imp2un
        self.dih0_imp3un = dih0_imp3un
        self.dih0_imp4un = dih0_imp4un
        self.dih0_native = dih0_native
        self.dih0_coef = dih0_coef
        self.dih0_type = dih0_type_str
        self.dih1_id = id_dih1
        self.dih1_iunit1 = dih1_iunit1
        self.dih1_iunit2 = dih1_iunit2
        self.dih1_imp1 = dih1_imp1
        self.dih1_imp2 = dih1_imp2
        self.dih1_imp3 = dih1_imp3
        self.dih1_imp4 = dih1_imp4
        self.dih1_imp1un = dih1_imp1un
        self.dih1_imp2un = dih1_imp2un
        self.dih1_imp3un = dih1_imp3un
        self.dih1_imp4un = dih1_imp4un
        self.dih1_native = dih1_native
        self.dih1_coef = dih1_coef
        self.dih1_type = dih1_type_str
        self.dih2_id = id_dih2
        self.dih2_iunit1 = dih2_iunit1
        self.dih2_iunit2 = dih2_iunit2
        self.dih2_imp1 = dih2_imp1
        self.dih2_imp2 = dih2_imp2
        self.dih2_imp3 = dih2_imp3
        self.dih2_imp4 = dih2_imp4
        self.dih2_imp1un = dih2_imp1un
        self.dih2_imp2un = dih2_imp2un
        self.dih2_imp3un = dih2_imp3un
        self.dih2_imp4un = dih2_imp4un
        self.dih2_native = dih2_native
        self.dih2_coef = dih2_coef
        self.dih2_type = dih2_type_str
        self.sectert = sectert
        self.nHB = nHB
        self.atoms1 = atoms1
        self.atoms2 = atoms2

class TertiaryStackDT(Ninfo):
    _counter = 0
    _counter_ang = 0
    _counter_dih = 0
    def __init__(self,id=None,id_ang1=None,id_ang2=None,id_dih0=None,id_dih1=None,id_dih2=None,
                 iunit1=None,iunit2=None, imp1=None,imp2=None,imp1un=None,imp2un=None,
                 native=None,factor=None,correct_mgo=None,coef=None,type_str=None,
                 ang1_iunit1=None,ang1_iunit2=None,
                 ang1_imp1=None,ang1_imp2=None,ang1_imp1un=None,ang1_imp2un=None,ang1_imp3=None,ang1_imp3un=None,
                 ang1_native=None,ang1_coef=None,ang1_type_str=None,
                 ang2_iunit1=None,ang2_iunit2=None,
                 ang2_imp1=None,ang2_imp2=None,ang2_imp1un=None,ang2_imp2un=None,ang2_imp3=None,ang2_imp3un=None,
                 ang2_native=None,ang2_coef=None,ang2_type_str=None,
                 dih0_imp1=None, dih0_imp2=None, dih0_imp3=None, dih0_imp4=None,dih0_iunit1=None,dih0_iunit2=None,
                 dih0_imp1un=None, dih0_imp2un=None, dih0_imp3un=None, dih0_imp4un=None,
                 dih0_native=None,dih0_coef=None,dih0_type_str=None,
                 dih1_imp1=None, dih1_imp2=None, dih1_imp3=None, dih1_imp4=None,dih1_iunit1=None,dih1_iunit2=None,
                 dih1_imp1un=None, dih1_imp2un=None, dih1_imp3un=None, dih1_imp4un=None,
                 dih1_native=None,dih1_coef=None,dih1_type_str=None,
                 dih2_imp1=None, dih2_imp2=None, dih2_imp3=None, dih2_imp4=None,dih2_iunit1=None,dih2_iunit2=None,
                 dih2_imp1un=None, dih2_imp2un=None, dih2_imp3un=None, dih2_imp4un=None,
                 dih2_native=None,dih2_coef=None,dih2_type_str=None, 
                 excess1=None, excess2=None):
        if id is None:
            HBondDT._counter += 1
            id = HBondDT._counter
        if id_ang1 is None:
            HBondDT._counter_ang += 1
            id_ang1 = HBondDT._counter_ang
        if id_ang2 is None:
            HBondDT._counter_ang += 1
            id_ang2 = HBondDT._counter_ang
        if id_dih0 is None:
            HBondDT._counter_dih += 1
            id_dih0 = HBondDT._counter_dih
        if id_dih1 is None:
            HBondDT._counter_dih += 1
            id_dih1 = HBondDT._counter_dih
        if id_dih2 is None:
            HBondDT._counter_dih += 1
            id_dih2 = HBondDT._counter_dih
        Ninfo.__init__(self,id=id,iunit1=iunit1,iunit2=iunit2, imp1=imp1,imp2=imp2,imp1un=imp1un,imp2un=imp2un,
                       native=native,factor=factor,correct_mgo=correct_mgo,coef=coef,type_str=type_str)
        self.ang1_id = id_ang1
        self.ang1_iunit1 = ang1_iunit1
        self.ang1_iunit2 = ang1_iunit2
        self.ang1_imp1 = ang1_imp1
        self.ang1_imp2 = ang1_imp2
        self.ang1_imp3 = ang1_imp3
        self.ang1_imp1un = ang1_imp1un
        self.ang1_imp2un = ang1_imp2un
        self.ang1_imp3un = ang1_imp3un
        self.ang1_native = ang1_native
        self.ang1_coef = ang1_coef
        self.ang1_type = ang1_type_str
        self.ang2_id = id_ang2
        self.ang2_iunit1 = ang2_iunit1
        self.ang2_iunit2 = ang2_iunit2
        self.ang2_imp1 = ang2_imp1
        self.ang2_imp2 = ang2_imp2
        self.ang2_imp3 = ang2_imp3
        self.ang2_imp1un = ang2_imp1un
        self.ang2_imp2un = ang2_imp2un
        self.ang2_imp3un = ang2_imp3un
        self.ang2_native = ang2_native
        self.ang2_coef = ang2_coef
        self.ang2_type = ang2_type_str
        self.dih0_id = id_dih0
        self.dih0_iunit1 = dih0_iunit1
        self.dih0_iunit2 = dih0_iunit2
        self.dih0_imp1 = dih0_imp1
        self.dih0_imp2 = dih0_imp2
        self.dih0_imp3 = dih0_imp3
        self.dih0_imp4 = dih0_imp4
        self.dih0_imp1un = dih0_imp1un
        self.dih0_imp2un = dih0_imp2un
        self.dih0_imp3un = dih0_imp3un
        self.dih0_imp4un = dih0_imp4un
        self.dih0_native = dih0_native
        self.dih0_coef = dih0_coef
        self.dih0_type = dih0_type_str
        self.dih1_id = id_dih1
        self.dih1_iunit1 = dih1_iunit1
        self.dih1_iunit2 = dih1_iunit2
        self.dih1_imp1 = dih1_imp1
        self.dih1_imp2 = dih1_imp2
        self.dih1_imp3 = dih1_imp3
        self.dih1_imp4 = dih1_imp4
        self.dih1_imp1un = dih1_imp1un
        self.dih1_imp2un = dih1_imp2un
        self.dih1_imp3un = dih1_imp3un
        self.dih1_imp4un = dih1_imp4un
        self.dih1_native = dih1_native
        self.dih1_coef = dih1_coef
        self.dih1_type = dih1_type_str
        self.dih2_id = id_dih2
        self.dih2_iunit1 = dih2_iunit1
        self.dih2_iunit2 = dih2_iunit2
        self.dih2_imp1 = dih2_imp1
        self.dih2_imp2 = dih2_imp2
        self.dih2_imp3 = dih2_imp3
        self.dih2_imp4 = dih2_imp4
        self.dih2_imp1un = dih2_imp1un
        self.dih2_imp2un = dih2_imp2un
        self.dih2_imp3un = dih2_imp3un
        self.dih2_imp4un = dih2_imp4un
        self.dih2_native = dih2_native
        self.dih2_coef = dih2_coef
        self.dih2_type = dih2_type_str
        self.excess1 = excess1
        self.excess2 = excess2
    
    
class NinfoSet(object): 
    def __init__(self) :
        self.bondlengths = []
        self.fenes = []
        self.bondangles = []
        self.dihedrals = []
        self.contacts = []
        self.LJs = []
        self.basestacks = []
        self.basepairs = []
        self.basestackDTs = []
        self.hbondDTs = []
        self.tertiarystackDTs = []
        self.max_unit = 0
        self.max_mp = 0
    
    def update_info(self):
        n = 0
        m = 0
        for bl in self.bondlengths :
            if bl.iunit1 > n :
                n = bl.iunit1
            if bl.iunit2 > n :
                n = bl.iunit2
            if bl.imp1 > m:
                m = bl.imp1
            if bl.imp2 > m:
                m = bl.imp2
        for fene in self.fenes :
            if fene.iunit1 > n :
                n = fene.iunit1
            if fene.iunit2 > n :
                n = fene.iunit2
            if fene.imp1 > m:
                m = fene.imp1
            if fene.imp2 > m:
                m = fene.imp2
        for ba in self.bondangles :
            if ba.iunit1 > n :
                n = ba.iunit1
            if ba.iunit2 > n :
                n = ba.iunit2
            if ba.imp1 > m:
                m = ba.imp1
            if ba.imp2 > m:
                m = ba.imp2
            if ba.imp3 > m:
                m = ba.imp3
        for dih in self.dihedrals :
            if dih.iunit1 > n :
                n = dih.iunit1 
            if dih.iunit2 > n :
                n = dih.iunit2
            if dih.imp1 > m:
                m = dih.imp1
            if dih.imp2 > m:
                m = dih.imp2
            if dih.imp3 > m:
                m = dih.imp3
            if dih.imp4 > m:
                m = dih.imp4
        for con in self.contacts :
            if con.iunit1 > n :
                n = con.iunit1
            if con.iunit2 > n :
                n = con.iunit2
            if con.imp1 > m:
                m = con.imp1
            if con.imp2 > m:
                m = con.imp2
        for LJ in self.LJs:
            if LJ.iunit1 > n:
                n = LJ.iunit1
            if LJ.iunit2 > n:
                n = LJ.iunit2
            if LJ.imp1 > m:
                m = LJ.imp1
            if LJ.imp2 > m:
                m = LJ.imp2
        for bs in self.basestacks :
            if bs.iunit1 > n :
                n = bs.iunit1 
            if bs.iunit2 > n :
                n = bs.iunit2
            if bs.imp1 > m:
                m = bs.imp1
            if bs.imp2 > m:
                m = bs.imp2
        for bp in self.basepairs :
            if bp.iunit1 > n :
                n = bp.iunit1 
            if bp.iunit2 > n :
                n = bp.iunit2
            if bp.imp1 > m:
                m = bp.imp1
            if bp.imp2 > m:
                m = bp.imp2
        for bs in self.basestackDTs :
            if bs.iunit1 > n :
                n = bs.iunit1 
            if bs.iunit2 > n :
                n = bs.iunit2
            if bs.imp1 > m:
                m = bs.imp1
            if bs.imp2 > m:
                m = bs.imp2
            if bs.dih1_imp1 > m:
                m = bs.dih1_imp1
            if bs.dih1_imp2 > m:
                m = bs.dih1_imp2
            if bs.dih1_imp3 > m:
                m = bs.dih1_imp3
            if bs.dih1_imp4 > m:
                m = bs.dih1_imp4
            if bs.dih2_imp1 > m:
                m = bs.dih2_imp1
            if bs.dih2_imp2 > m:
                m = bs.dih2_imp2
            if bs.dih2_imp3 > m:
                m = bs.dih2_imp3
            if bs.dih2_imp4 > m:
                m = bs.dih2_imp4
        for hb in self.hbondDTs :
            if hb.iunit1 > n :
                n = hb.iunit1 
            if hb.iunit2 > n :
                n = hb.iunit2
            if hb.imp1 > m:
                m = hb.imp1
            if hb.imp2 > m:
                m = hb.imp2
            if hb.ang1_imp1 > m:
                m = hb.ang1_imp1
            if hb.ang1_imp2 > m:
                m = hb.ang1_imp2
            if hb.ang1_imp3 > m:
                m = hb.ang1_imp3
            if hb.ang2_imp1 > m:
                m = hb.ang2_imp1
            if hb.ang2_imp2 > m:
                m = hb.ang2_imp2
            if hb.ang2_imp3 > m:
                m = hb.ang2_imp3
            if hb.dih0_imp1 > m:
                m = hb.dih0_imp1
            if hb.dih0_imp2 > m:
                m = hb.dih0_imp2
            if hb.dih0_imp3 > m:
                m = hb.dih0_imp3
            if hb.dih0_imp4 > m:
                m = hb.dih0_imp4
            if hb.dih1_imp1 > m:
                m = hb.dih1_imp1
            if hb.dih1_imp2 > m:
                m = hb.dih1_imp2
            if hb.dih1_imp3 > m:
                m = hb.dih1_imp3
            if hb.dih1_imp4 > m:
                m = hb.dih1_imp4
            if hb.dih2_imp1 > m:
                m = hb.dih2_imp1
            if hb.dih2_imp2 > m:
                m = hb.dih2_imp2
            if hb.dih2_imp3 > m:
                m = hb.dih2_imp3
            if hb.dih2_imp4 > m:
                m = hb.dih2_imp4
        for tst in self.tertiarystackDTs :
            if tst.iunit1 > n :
                n = tst.iunit1 
            if tst.iunit2 > n :
                n = tst.iunit2
            if tst.imp1 > m:
                m = tst.imp1
            if tst.imp2 > m:
                m = tst.imp2
            if tst.ang1_imp1 > m:
                m = tst.ang1_imp1
            if tst.ang1_imp2 > m:
                m = tst.ang1_imp2
            if tst.ang1_imp3 > m:
                m = tst.ang1_imp3
            if tst.ang2_imp1 > m:
                m = tst.ang2_imp1
            if tst.ang2_imp2 > m:
                m = tst.ang2_imp2
            if tst.ang2_imp3 > m:
                m = tst.ang2_imp3
            if tst.dih0_imp1 > m:
                m = tst.dih0_imp1
            if tst.dih0_imp2 > m:
                m = tst.dih0_imp2
            if tst.dih0_imp3 > m:
                m = tst.dih0_imp3
            if tst.dih0_imp4 > m:
                m = tst.dih0_imp4
            if tst.dih1_imp1 > m:
                m = tst.dih1_imp1
            if tst.dih1_imp2 > m:
                m = tst.dih1_imp2
            if tst.dih1_imp3 > m:
                m = tst.dih1_imp3
            if tst.dih1_imp4 > m:
                m = tst.dih1_imp4
            if tst.dih2_imp1 > m:
                m = tst.dih2_imp1
            if tst.dih2_imp2 > m:
                m = tst.dih2_imp2
            if tst.dih2_imp3 > m:
                m = tst.dih2_imp3
            if tst.dih2_imp4 > m:
                m = tst.dih2_imp4
        self.max_unit = n
        self.max_mp = m
        
    def dict_of_ninfoset_by_unit(self):
        di = {}
        self.update_info()
        for i in range(1, self.max_unit+1):
            for j in range(i, self.max_unit+1):
                ns = NinfoSet()
                di[(i,j)] = ns
        for bl in self.bondlengths:
            di[(bl.iunit1,bl.iunit1)].bondlengths.append(bl)
        for fene in self.fenes:
            di[(fene.iunit1,fene.iunit1)].fenes.append(fene)
        for ba in self.bondangles:
            di[(ba.iunit1,ba.iunit1)].bondangles.append(ba)
        for dih in self.dihedrals:
            di[(dih.iunit1,dih.iunit1)].dihedrals.append(dih)
        for con in self.contacts:
            di[(con.iunit1,con.iunit2)].contacts.append(con)
        for LJ in self.LJs:
            di[(LJ.iunit1,LJ.iunit2)].LJs.append(LJ)
        for bp in self.basepairs:
            di[(bp.iunit1,bp.iunit2)].basepairs.append(bp)
        for bs in self.basestacks:
            di[(bs.iunit1,bs.iunit2)].basestacks.append(bs)
        for bs in self.basestackDTs:
            di[(bs.iunit1,bs.iunit2)].basestackDTs.append(bs)
        for hb in self.hbondDTs:
            di[(hb.iunit1,hb.iunit2)].hbondDTs.append(hb)
        for tst in self.tertiarystackDTs:
            di[(tst.iunit1,tst.iunit2)].tertiarystackDTs.append(tst)
            
        for n in di:
            di[n].update_info()
        return di
        
        
    def get_contacts_by_unit(self,i,j):
        subset = []
        for con in self.contacts :
            if con.iunit1 == i and con.iunit2 == j :
                subset.append(con)
            elif con.iunit1 == j and con.iunit2 == i :
                subset.append(con)
        return subset
        
    def get_LJs_by_unit(self,i,j):
        subset = []
        for LJ in self.LJs :
            if LJ.iunit1 == i and LJ.iunit2 == j :
                subset.append(LJ)
            elif LJ.iunit1 == j and LJ.iunit2 == i :
                subset.append(LJ)
        return subset
    
    def get_basepairs_by_unit(self,i,j):
        subset = []
        for con in self.basepairs :
            if con.iunit1 == i and con.iunit2 == j :
                subset.append(con)
            elif con.iunit1 == j and con.iunit2 == i :
                subset.append(con)
        return subset
    
    def reassign_id(self):
        for newid, bl in enumerate(self.bondlengths) :
            bl.id = newid + 1
        for newid, fene in enumerate(self.fenes) :
            fene.id = newid + 1
        for newid, ba in enumerate(self.bondangles) :
            ba.id = newid + 1
        for newid, con in enumerate(self.contacts) :
            con.id = newid + 1
        for newid, LJ in enumerate(self.LJs) :
            LJ.id = newid + 1
        for newid, bp in enumerate(self.basepairs) :
            bp.id = newid + 1
        for newid, bs in enumerate(self.basestacks) :
            bs.id = newid + 1
        ## Need to update angle and dihedral associated IDs
        #for newid, bs in enumerate(self.basestackDTs) :
        #    bs.id = newid + 1
        for newid, hb in enumerate(self.hbondDTs):
            hb.id = newid + 1
            hb.ang1_id = 2*(newid+1) - 1
            hb.ang2_id = 2*(newid+1)
            hb.dih0_id = 3*(newid+1) - 2
            hb.dih1_id = 3*(newid+1) - 1
            hb.dih2_id = 3*(newid+1)
        for newid, tst in enumerate(self.tertiarystackDTs):
            tst.id = newid + 1
            tst.ang1_id = 2*(newid+1) - 1
            tst.ang2_id = 2*(newid+1)
            tst.dih0_id = 3*(newid+1) - 2
            tst.dih1_id = 3*(newid+1) - 1
            tst.dih2_id = 3*(newid+1)
            
    def sort_by_mp(self):
        # hbond
        newlist = []
        for imp1 in range(1,self.max_mp+1):
            for imp2 in range(imp1,self.max_mp+1):
                for hb in self.hbondDTs:
                    if ((hb.imp1 == imp1 and hb.imp2 == imp2) or (hb.imp1 == imp2 and hb.imp2 == imp1)):
                        newlist.append( copy.deepcopy( hb ))
        self.hbondDTs = newlist

        # Tertiarystack
        newlist = []
        for imp1 in range(1,self.max_mp+1):
            for imp2 in range(imp1,self.max_mp+1):
                for tst in self.tertiarystackDTs:
                    if ((tst.imp1 == imp1 and tst.imp2 == imp2) or (tst.imp1 == imp2 and tst.imp2 == imp1)):
                        newlist.append( copy.deepcopy( tst ))
        self.tertiarystackDTs = newlist



