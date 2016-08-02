#!/usr/bin/env python
'''
@author: Naoto Hori
'''

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
        
#**      ibd iunit1-iunit2   imp1 - imp2 imp1un-imp2un      bd_nat    factor_bd  correct_mgo      coef_bd
#bond      1      1      1      1      2      1      2       3.8531       1.0000       1.0000     110.4000 pp
class BondLength(Ninfo):
    _counter = 0
    def __init__(self,id=None,iunit1=None,iunit2=None, imp1=None,imp2=None,imp1un=None,imp2un=None,
                 native=None,factor=None,correct_mgo=None,coef=None,type_str=None):
        if id == None:
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
        if id == None:
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
        if id == None:
            BondAngle._counter += 1
            id = BondAngle._counter
        Ninfo.__init__(self,id=id,iunit1=iunit1,iunit2=iunit2, imp1=imp1,imp2=imp2,imp1un=imp1un,imp2un=imp2un,
                       native=native,factor=factor,correct_mgo=correct_mgo,coef=coef,type_str=type_str)
        self.imp3 = imp3
        self.imp3un = imp3un
    
    def show(self):
        Ninfo.show(self)
        print 'imp3',self.imp3
        print 'imp3un',self.imp3un

class Dihedral(Ninfo):
#**     idih iunit1-iunit2   imp1 - imp2 - imp3 - imp4 imp1un-imp2un-imp3un-imp4un      dih_nat   factor_dih  correct_mgo   coef_dih_1   coef_dih_3
#dihd      1      1      1      1      2      3      4      1      2      3      4    -132.8457       0.0000       1.0000       0.0000       0.0000 pppp
    _counter = 0
    def __init__(self,id=None,iunit1=None,iunit2=None, imp1=None,imp2=None,imp1un=None,imp2un=None,
                 imp3=None,imp3un=None,imp4=None,imp4un=None,coef_3=None,
                 native=None,factor=None,correct_mgo=None,coef=None,type_str=None):
        if id == None:
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
        print 'imp3',self.imp3
        print 'imp4',self.imp4
        print 'imp3un',self.imp3un
        print 'imp4un',self.imp4un
        print 'coef_3',self.coef_3

class Contact(Ninfo):
    _counter = 0
    def __init__(self,id=None,iunit1=None,iunit2=None, imp1=None,imp2=None,imp1un=None,imp2un=None,dummy=None,
                 native=None,factor=None,correct_mgo=None,coef=None,type_str=None):
        if id == None:
            Contact._counter += 1
            id = Contact._counter
        Ninfo.__init__(self,id=id,iunit1=iunit1,iunit2=iunit2, imp1=imp1,imp2=imp2,imp1un=imp1un,imp2un=imp2un,
                       native=native,factor=factor,correct_mgo=correct_mgo,coef=coef,type_str=type_str)
        self.dummy = dummy

class LJ(Ninfo):
    _counter = 0 
    def __init__(self,id=None,iunit1=None,iunit2=None, imp1=None,imp2=None,imp1un=None,imp2un=None,
                 native=None,factor=None,correct_mgo=None,coef=None,type_str=None):
        if id == None:
            LJ._counter += 1
            id = LJ._counter
        Ninfo.__init__(self,id=id,iunit1=iunit1,iunit2=iunit2, imp1=imp1,imp2=imp2,imp1un=imp1un,imp2un=imp2un,
                       native=native,factor=factor,correct_mgo=correct_mgo,coef=coef,type_str=type_str)
    
class BasePair(Ninfo):
    _counter = 0
    def __init__(self,id=None,iunit1=None,iunit2=None, imp1=None,imp2=None,imp1un=None,imp2un=None,dummy=None,nhb=None,
                 native=None,factor=None,correct_mgo=None,coef=None,type_str=None):
        if id == None:
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
        if id == None:
            BaseStack._counter += 1
            id = BaseStack._counter
        Ninfo.__init__(self,id=id,iunit1=iunit1,iunit2=iunit2, imp1=imp1,imp2=imp2,imp1un=imp1un,imp2un=imp2un,
                       native=native,factor=factor,correct_mgo=correct_mgo,coef=coef,type_str=type_str)
        self.dummy = dummy

class BaseStackDT13(Ninfo):
    _counter = 0
    _counter_dih = 0
    def __init__(self,id=None,id_dih1=None,id_dih2=None,
                 iunit1=None,iunit2=None, imp1=None,imp2=None,imp1un=None,imp2un=None,
                 native=None,factor=None,correct_mgo=None,coef=None,type_str=None,
                 dih1_imp1=None, dih1_imp2=None, dih1_imp3=None, dih1_imp4=None,dih1_iunit1=None,dih1_iunit2=None,
                 dih1_imp1un=None, dih1_imp2un=None, dih1_imp3un=None, dih1_imp4un=None,
                 dih1_native=None,dih1_coef=None,dih1_type_str=None,
                 dih2_imp1=None, dih2_imp2=None, dih2_imp3=None, dih2_imp4=None,dih2_iunit1=None,dih2_iunit2=None,
                 dih2_imp1un=None, dih2_imp2un=None, dih2_imp3un=None, dih2_imp4un=None,
                 dih2_native=None,dih2_coef=None,dih2_type_str=None):
        if id == None:
            BaseStackDT13._counter += 1
            id = BaseStackDT13._counter
        if id_dih1 == None:
            BaseStackDT13._counter_dih += 1
            id_dih1 = BaseStackDT13._counter_dih
        if id_dih2 == None:
            BaseStackDT13._counter_dih += 1
            id_dih2 = BaseStackDT13._counter_dih
        Ninfo.__init__(self,id=id,iunit1=iunit1,iunit2=iunit2, imp1=imp1,imp2=imp2,imp1un=imp1un,imp2un=imp2un,
                       native=native,factor=factor,correct_mgo=correct_mgo,coef=coef,type_str=type_str)
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
    
class NinfoSet(object): 
    def __init__(self) :
        self.bondlengths = []
        self.fenes = []
        self.bondangles = []
        self.dihedrals = []
        self.contacts = []
        self.LJs = []
        self.basestacks = []
        self.basestackDT13s = []
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
        for fene in self.fenes :
            if fene.iunit1 > n :
                n = fene.iunit1
            if fene.iunit2 > n :
                n = fene.iunit2
        for ba in self.bondangles :
            if ba.iunit1 > n :
                n = ba.iunit1
            if ba.iunit2 > n :
                n = ba.iunit2
        for dih in self.dihedrals :
            if dih.iunit1 > n :
                n = dih.iunit1 
            if dih.iunit2 > n :
                n = dih.iunit2
        for con in self.contacts :
            if con.iunit1 > n :
                n = con.iunit1
            if con.iunit2 > n :
                n = con.iunit2
        for LJ in self.LJs:
            if LJ.iunit1 > n:
                n = LJ.iunit1
            if LJ.iunit2 > n:
                n = LJ.iunit2
        for bs in self.basestacks :
            if bs.iunit1 > n :
                n = bs.iunit1 
            if bs.iunit2 > n :
                n = bs.iunit2
        for bs in self.basestackDT13s :
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
        for bs in self.basestackDT13s:
            di[(bs.iunit1,bs.iunit2)].basestackDT13s.append(bs)
            
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
        for newid, bs in enumerate(self.basestackDT13s) :
            bs.id = newid + 1
            
