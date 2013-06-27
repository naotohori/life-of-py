#!/usr/bin/env python
'''
@author: Naoto Hori
'''

class Ninfo(object):
    def __init__(self):
        self.id = None
        self.iunit1 = None
        self.iunit2 = None
        self.imp1 = None
        self.imp2 = None
        self.imp1un = None
        self.imp2un = None
        self.native = None
        self.factor = None
        self.correct_mgo = None
        self.coef = None
        self.type = None

    def show(self):
        print 'id',self.id
        print 'iunit1',self.iunit1
        print 'iunit2',self.iunit2
        print 'imp1',self.imp1
        print 'imp2',self.imp2
        print 'imp1un',self.imp1un
        print 'imp2un',self.imp2un
        print 'native',self.native
        print 'factor',self.factor
        print 'correct_mgo',self.correct_mgo
        print 'coef',self.coef
        print 'type',self.type
        
class BondLength(Ninfo):
#**      ibd iunit1-iunit2   imp1 - imp2 imp1un-imp2un      bd_nat    factor_bd  correct_mgo      coef_bd
#bond      1      1      1      1      2      1      2       3.8531       1.0000       1.0000     110.4000 pp
    pass
    
class BondAngle(Ninfo):
#**      iba iunit1-iunit2   imp1 - imp2 - imp3 imp1un-imp2un-imp3un      ba_nat    factor_ba  correct_mgo      coef_ba
#angl      1      1      1      1      2      3      1      2      3     114.0257       0.0000      1.0000       0.0000 ppp
    imp3 = None
    imp3un = None
    
    def show(self):
        Ninfo.show(self)
        print 'imp3',self.imp3
        print 'imp3un',self.imp3un

class Aicg13(Ninfo):
#**      iba iunit1-iunit2   imp1 - imp2 - imp3 imp1un-imp2un-imp3un  aicg13_nat  factor_aicg13  correct_mgo  coef_aicg13_gauss wid_aicg13_gauss
#aicg13      1      1      1      1      2      3      1      2      3       6.4310       1.0000       1.0000       1.1594       0.1500 ppp
    imp3 = None
    imp3un = None
    wid = None
    
    def show(self):
        Ninfo.show(self)
        print 'imp3',self.imp3
        print 'imp3un',self.imp3un
        print 'wid',self.wid

class Dihedral(Ninfo):
#**     idih iunit1-iunit2   imp1 - imp2 - imp3 - imp4 imp1un-imp2un-imp3un-imp4un      dih_nat   factor_dih  correct_mgo   coef_dih_1   coef_dih_3
#dihd      1      1      1      1      2      3      4      1      2      3      4    -132.8457       0.0000       1.0000       0.0000       0.0000 pppp
    imp3 = None
    imp4 = None
    imp3un = None
    imp4un = None
    coef_3 = None
    
    def show(self):
        Ninfo.show(self)
        print 'imp3',self.imp3
        print 'imp4',self.imp4
        print 'imp3un',self.imp3un
        print 'imp4un',self.imp4un
        print 'coef_3',self.coef_3

class AicgDih(Ninfo):
#**     idih iunit1-iunit2   imp1 - imp2 - imp3 - imp4 imp1un-imp2un-imp3un-imp4un   dih_nat factor_aicg14  correct_mgo  coef_dih_gauss  wid_dih_gauss
#aicgdih      1      1      1      1      2      3      4      1      2      3      4    -132.8457       1.0000       1.0000       0.3013       0.1500 pppp
    imp3 = None
    imp4 = None
    imp3un = None
    imp4un = None
    wid = None
    
    def show(self):
        Ninfo.show(self)
        print 'imp3',self.imp3
        print 'imp4',self.imp4
        print 'imp3un',self.imp3un
        print 'imp4un',self.imp4un
        print 'wid',self.wid

class Contact(Ninfo):
    dummy = None
    
class BasePair(Ninfo):
    dummy = None
    nhb = None
    
class BaseStack(Ninfo):
    dummy = None
    
class NinfoSet(object): 
    def __init__(self) :
        self.bondlengths = []
        self.bondangles = []
        self.aicg13s = []
        self.dihedrals = []
        self.aicgdihs = []
        self.contacts = []
        self.basestacks = []
        self.basepairs = []
        self.max_unit = 0
    
    def update_info(self):
        n = 0
#        n = self.max_unit
        for bl in self.bondlengths :
            if bl.iunit1 > n :
                n = bl.iunit1
            if bl.iunit2 > n :
                n = bl.iunit2
        for ba in self.bondangles :
            if ba.iunit1 > n :
                n = ba.iunit1
            if ba.iunit2 > n :
                n = ba.iunit2
        for a13 in self.aicg13s :
            if a13.iunit1 > n :
                n = a13.iunit1
            if a13.iunit2 > n :
                n = a13.iunit2
        for dih in self.dihedrals :
            if dih.iunit1 > n :
                n = dih.iunit1 
            if dih.iunit2 > n :
                n = dih.iunit2
        for adih in self.aicgdihs :
            if adih.iunit1 > n :
                n = adih.iunit1 
            if adih.iunit2 > n :
                n = adih.iunit2
        for con in self.contacts :
            if con.iunit1 > n :
                n = con.iunit1
            if con.iunit2 > n :
                n = con.iunit2
        for bs in self.basestacks :
            if bs.iunit1 > n :
                n = bs.iunit1 
            if bs.iunit2 > n :
                n = bs.iunit2
        for bp in self.basepairs :
            if bp.iunit1 > n :
                n = bp.iunit1 
            if bp.iunit2 > n :
                n = bp.iunit2
        self.max_unit = n
        
    def dict_of_ninfoset_by_unit(self):
        di = {}
        self.update_info()
        for i in xrange(1, self.max_unit+1):
            for j in xrange(i, self.max_unit+1):
                ns = NinfoSet()
                di[(i,j)] = ns
        for bl in self.bondlengths:
            di[(bl.iunit1,bl.iunit1)].bondlengths.append(bl)
        for ba in self.bondangles:
            di[(ba.iunit1,ba.iunit1)].bondangles.append(ba)
        for a13 in self.aicg13s:
            di[(a13.iunit1,a13.iunit1)].aicg13s.append(a13)
        for dih in self.dihedrals:
            di[(dih.iunit1,dih.iunit1)].dihedrals.append(dih)
        for adih in self.aicgdihs:
            di[(adih.iunit1,adih.iunit1)].aicgdihs.append(adih)
        for con in self.contacts:
            di[(con.iunit1,con.iunit2)].contacts.append(con)
        for bp in self.basepairs:
            di[(bp.iunit1,bp.iunit2)].basepairs.append(bp)
        for bs in self.basestacks:
            di[(bs.iunit1,bs.iunit2)].basestacks.append(bs)
            
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
        for newid, ba in enumerate(self.bondangles) :
            ba.id = newid + 1
        for newid, a13 in enumerate(self.aicg13s) :
            a13.id = newid + 1
        for newid, adih in enumerate(self.aicgdihs) :
            adih.id = newid + 1
        for newid, con in enumerate(self.contacts) :
            con.id = newid + 1
        for newid, bp in enumerate(self.basepairs) :
            bp.id = newid + 1
        for newid, bs in enumerate(self.basestacks) :
            bs.id = newid + 1
            