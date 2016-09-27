#!/usr/bin/env python

seq = 'UAUUUGUUUGCUUUUUUAUGCUUUUUUUUUUUUUUUUUUU'
n_nt = len(seq)
#n_mp = 3 * n_nt - 1

from cafysis.para.rnaAform import ARNA
from cafysis.para.rnaDT13 import DT13
from cafysis.elements.ninfo import NinfoSet, BondLength, BondAngle, BaseStackDT

ns = NinfoSet()

##############
## Bond     ##
##############

def bl_native(i):
    if seq[i-1] == 'A':
        native = ARNA.BL_SA
        type_str = 'SA'
    elif seq[i-1] == 'U':
        native = ARNA.BL_SU
        type_str = 'SU'
    elif seq[i-1] == 'G':
        native = ARNA.BL_SG
        type_str = 'SG'
    elif seq[i-1] == 'C':
        native = ARNA.BL_SC
        type_str = 'SC'
    return native, type_str

# for the first nt
# S - B
native, type_str = bl_native(1)
bl = BondLength(iunit1=1,iunit2=1,imp1=1,imp2=2,imp1un=1,imp2un=2,
                native=native,factor=1.0,correct_mgo=1.0,coef=DT13.BL_SB,type_str=type_str)
ns.bondlengths.append(bl)

# for the second through last nt
for i in range(2,n_nt+1):
    imp_P = 3 * (i-1)
    imp_S0 = imp_P - 2
    imp_S = imp_P + 1
    imp_B = imp_S + 1
    # S0 - P (from previous nt)
    bl = BondLength(iunit1=1,iunit2=1,imp1=imp_S0,imp2=imp_P,imp1un=imp_S0,imp2un=imp_P,
                    native=ARNA.BL_SP,factor=1.0,correct_mgo=1.0,coef=DT13.BL_SP,type_str='SP')
    ns.bondlengths.append(bl)

    # P - S
    bl = BondLength(iunit1=1,iunit2=1,imp1=imp_P,imp2=imp_S,imp1un=imp_P,imp2un=imp_S,
                    native=ARNA.BL_PS,factor=1.0,correct_mgo=1.0,coef=DT13.BL_PS,type_str='PS')
    ns.bondlengths.append(bl)

    # S - B
    native, type_str = bl_native(i)
    bl = BondLength(iunit1=1,iunit2=1,imp1=imp_S,imp2=imp_B,imp1un=imp_S,imp2un=imp_B,
                    native=native,factor=1.0,correct_mgo=1.0,coef=DT13.BL_SB,type_str=type_str)
    ns.bondlengths.append(bl)


##############
## Angle    ##
##############

def ba_BSP_native(i):
    if seq[i-1] == 'A':
        native = ARNA.BA_ASP
        type_str = 'ASP'
    elif seq[i-1] == 'U':
        native = ARNA.BA_USP
        type_str = 'USP'
    elif seq[i-1] == 'G':
        native = ARNA.BA_GSP
        type_str = 'GSP'
    elif seq[i-1] == 'C':
        native = ARNA.BA_CSP
        type_str = 'CSP'
    return native, type_str
def ba_PSB_native(i):
    if seq[i-1] == 'A':
        native = ARNA.BA_PSA
        type_str = 'PSA'
    elif seq[i-1] == 'U':
        native = ARNA.BA_PSU
        type_str = 'PSU'
    elif seq[i-1] == 'G':
        native = ARNA.BA_PSG
        type_str = 'PSG'
    elif seq[i-1] == 'C':
        native = ARNA.BA_PSC
        type_str = 'PSC'
    return native, type_str

# for the second nt
# B0 - S0 - P
native, type_str = ba_BSP_native(1)
ba = BondAngle(iunit1=1,iunit2=1,imp1=2,imp2=1,imp3=3,imp1un=2,imp2un=1,imp3un=3,
                native=native,factor=1.0,correct_mgo=1.0,coef=DT13.BA_BSP,type_str=type_str)
ns.bondangles.append(ba)
# S0 - P - S
ba = BondAngle(iunit1=1,iunit2=1,imp1=1,imp2=3,imp3=4,imp1un=1,imp2un=3,imp3un=4,
                native=ARNA.BA_SPS,factor=1.0,correct_mgo=1.0,coef=DT13.BA_SPS,type_str='SPS')
ns.bondangles.append(ba)
# P - S - B
native, type_str = ba_PSB_native(2)
ba = BondAngle(iunit1=1,iunit2=1,imp1=3,imp2=4,imp3=5,imp1un=3,imp2un=4,imp3un=5,
                native=native,factor=1.0,correct_mgo=1.0,coef=DT13.BA_PSB,type_str=type_str)
ns.bondangles.append(ba)

# for the third through last nt
for i in range(3,n_nt+1):
    imp_P = 3 * (i-1)
    imp_P0 = imp_P - 3
    imp_S0 = imp_P - 2
    imp_B0 = imp_P - 1
    imp_S = imp_P + 1
    imp_B = imp_S + 1
    # P0 - S0 - P
    ba = BondAngle(iunit1=1,iunit2=1,imp1=imp_P0,imp2=imp_S0,imp3=imp_P,imp1un=imp_P0,imp2un=imp_S0,imp3un=imp_P,
                    native=ARNA.BA_PSP,factor=1.0,correct_mgo=1.0,coef=DT13.BA_PSP,type_str='PSP')
    ns.bondangles.append(ba)
    # B0 - S0 - P
    native, type_str = ba_BSP_native(i-1)
    ba = BondAngle(iunit1=1,iunit2=1,imp1=imp_B0,imp2=imp_S0,imp3=imp_P,imp1un=imp_B0,imp2un=imp_S0,imp3un=imp_P,
                    native=native,factor=1.0,correct_mgo=1.0,coef=DT13.BA_BSP,type_str=type_str)
    ns.bondangles.append(ba)
    # S0 - P - S
    ba = BondAngle(iunit1=1,iunit2=1,imp1=imp_S0,imp2=imp_P,imp3=imp_S,imp1un=imp_S0,imp2un=imp_P,imp3un=imp_S,
                    native=ARNA.BA_SPS,factor=1.0,correct_mgo=1.0,coef=DT13.BA_SPS,type_str='SPS')
    ns.bondangles.append(ba)
    # P - S - B
    native, type_str = ba_PSB_native(i)
    ba = BondAngle(iunit1=1,iunit2=1,imp1=imp_P,imp2=imp_S,imp3=imp_B,imp1un=imp_P,imp2un=imp_S,imp3un=imp_B,
                    native=native,factor=1.0,correct_mgo=1.0,coef=DT13.BA_PSB,type_str=type_str)
    ns.bondangles.append(ba)

##############
## Stack    ##
##############
def bs_native(i,j):
    if seq[i-1] == 'A':
        if seq[j-1] == 'A':
            native = ARNA.ST_AA
            type_str = 'A-A'
        elif seq[j-1] == 'U':
            native = ARNA.ST_AU
            type_str = 'A-U'
        elif seq[j-1] == 'G':
            native = ARNA.ST_AG
            type_str = 'A-G'
        elif seq[j-1] == 'C':
            native = ARNA.ST_AC
            type_str = 'A-C'
    elif seq[i-1] == 'U':
        if seq[j-1] == 'A':
            native = ARNA.ST_UA
            type_str = 'U-A'
        elif seq[j-1] == 'U':
            native = ARNA.ST_UU
            type_str = 'U-U'
        elif seq[j-1] == 'G':
            native = ARNA.ST_UG
            type_str = 'U-G'
        elif seq[j-1] == 'C':
            native = ARNA.ST_UC
            type_str = 'U-C'
    elif seq[i-1] == 'G':
        if seq[j-1] == 'A':
            native = ARNA.ST_GA
            type_str = 'G-A'
        elif seq[j-1] == 'U':
            native = ARNA.ST_GU
            type_str = 'G-U'
        elif seq[j-1] == 'G':
            native = ARNA.ST_GG
            type_str = 'G-G'
        elif seq[j-1] == 'C':
            native = ARNA.ST_GC
            type_str = 'G-C'
    elif seq[i-1] == 'C':
        if seq[j-1] == 'A':
            native = ARNA.ST_CA
            type_str = 'C-A'
        elif seq[j-1] == 'U':
            native = ARNA.ST_CU
            type_str = 'C-U'
        elif seq[j-1] == 'G':
            native = ARNA.ST_CG
            type_str = 'C-G'
        elif seq[j-1] == 'C':
            native = ARNA.ST_CC
            type_str = 'C-C'
    return native, type_str

# for the second through one before last nt
for i in range(3,n_nt):
    imp_P2 = 3 * (i-1)
    imp_P1 = imp_P2 - 3
    imp_S1 = imp_P2 - 2
    imp_B1 = imp_P2 - 1
    imp_S2 = imp_P2 + 1
    imp_B2 = imp_P2 + 2
    imp_P3 = imp_P2 + 3

    # dist
    native, type_str = bs_native(i-1,i)
    bs = BaseStackDT(iunit1=1,iunit2=1,imp1=imp_B1,imp2=imp_B2,imp1un=imp_B1,imp2un=imp_B2,
                       native=native, factor=0.0,correct_mgo=1.0,coef=DT13.ST_DIST,type_str=type_str,
                       dih1_imp1=imp_P1, dih1_imp2=imp_S1, dih1_imp3=imp_P2, dih1_imp4=imp_S2,dih1_iunit1=1,dih1_iunit2=1,
                       dih1_imp1un=imp_P1, dih1_imp2un=imp_S1, dih1_imp3un=imp_P2, dih1_imp4un=imp_S2,
                       dih1_native=ARNA.DIH_PSPS,dih1_coef=DT13.ST_DIH,dih1_type_str='PSPS',
                       dih2_imp1=imp_S1, dih2_imp2=imp_P2, dih2_imp3=imp_S2, dih2_imp4=imp_P3,dih2_iunit1=1,dih2_iunit2=1,
                       dih2_imp1un=imp_S1, dih2_imp2un=imp_P2, dih2_imp3un=imp_S2, dih2_imp4un=imp_P3,
                       dih2_native=ARNA.DIH_SPSP,dih2_coef=DT13.ST_DIH,dih2_type_str='SPSP')
    ns.basestackDTs.append(bs)

from cafysis.file_io.ninfo import NinfoFile
nf = NinfoFile('ninfo_ssRNA13.ninfo')
nf.open_to_write()
nf.write_all(ns)
nf.close()

