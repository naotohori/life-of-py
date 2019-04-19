#!/usr/bin/env python

import sys
import math

from cafysis.file_io.pdb import PdbFile
from cafysis.file_io.ninfo import NinfoFile
from cafysis.para.rnaAform import ARNA
from cafysis.para.rnaDT13 import DT13
from cafysis.elements.ninfo import NinfoSet, BondLength, BondAngle, BaseStackDT, HBondDT

if len(sys.argv) != 4:
    print('Usage: SCRIPT [cg pdb] [hb list file] [output ninfo]')
    sys.exit(2)

f_in = PdbFile(sys.argv[1])
f_in.open_to_read()
chains = f_in.read_all()
f_in.close()

#seq = 'UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU'
#n_nt = len(seq)
#n_mp = 3 * n_nt - 1

if len(chains) > 1:
    print('Error: currently more than one chain can not be processed.')
    sys.exit(2)

chain = chains[0]
n_nt = chain.num_res()
print('#nt: ', n_nt)

seq = []
for r in chain.residues:
    # "RA " ---> "A"
    seq.append(r.atoms[0].res_name.strip()[1])
print('Sequence:')
print(seq)

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

##############
## H-bond   ##
##############
hblist = []
for l in open(sys.argv[2],'r'):
    if l.find('#') != -1:
        continue
    lsp = l.split()

    ## Check
    if lsp[0] == 'CAN':
        if lsp[2] != 'B' or lsp[4] != 'B':
            print('Error: Canonical base pair should be by B and B')
            sys.exit(2)
    elif lsp[0] == 'NON':
        pass
    else:
        print('Error: unknown H-bond type')
        sys.exit(2)

    hblist.append((lsp[0],int(lsp[1]),lsp[2],int(lsp[3]),lsp[4],int(lsp[5])))

def hb_ARNA_native(i,j):
    if seq[i-1] == 'A' and seq[j-1] == 'U':
        dist = ARNA.HB_AU
        ang1 = ARNA.HBA_SAU
        ang2 = ARNA.HBA_SUA
        dih0 = ARNA.HBD_SAUS
        dih1 = ARNA.HBD_PSAU
        dih2 = ARNA.HBD_PSUA
        nHB = 2
    elif seq[i-1] == 'U' and seq[j-1] == 'A':
        dist = ARNA.HB_AU
        ang1 = ARNA.HBA_SUA
        ang2 = ARNA.HBA_SAU
        dih0 = ARNA.HBD_SAUS
        dih1 = ARNA.HBD_PSUA
        dih2 = ARNA.HBD_PSAU
        nHB = 2
    elif seq[i-1] == 'G' and seq[j-1] == 'C':
        dist = ARNA.HB_GC
        ang1 = ARNA.HBA_SGC
        ang2 = ARNA.HBA_SCG
        dih0 = ARNA.HBD_SGCS
        dih1 = ARNA.HBD_PSGC
        dih2 = ARNA.HBD_PSCG
        nHB = 3
    elif seq[i-1] == 'C' and seq[j-1] == 'G':
        dist = ARNA.HB_GC
        ang1 = ARNA.HBA_SCG
        ang2 = ARNA.HBA_SGC
        dih0 = ARNA.HBD_SGCS
        dih1 = ARNA.HBD_PSCG
        dih2 = ARNA.HBD_PSGC
        nHB = 3
    else:
        print('Canonical basepair should be A-U or G-C: ',i,j)
        sys.exit(2)
    return dist, ang1, ang2, dih0, dih1, dih2, nHB


for c in hblist:
    nt_1 = c[1]
    mp_1 = c[2]
    nt_2 = c[3]
    mp_2 = c[4]
    nHB  = c[5]

    # For both CAN and NON, 
    if mp_1 == 'S':
        imp_1 = 2 + 3 * (nt_1 - 1) - 1  # S
        imp_3 = 2 + 3 * (nt_1 - 1) + 1  # P
        imp_5 = 2 + 3 * (nt_1 - 1) + 2  # S
    elif mp_1 == 'B':
        imp_1 = 2 + 3 * (nt_1 - 1)      # B
        imp_3 = 2 + 3 * (nt_1 - 1) - 1  # S
        imp_5 = 2 + 3 * (nt_1 - 1) + 1  # P
    elif mp_1 == 'P':
        imp_1 = 2 + 3 * (nt_1 - 1) - 2  # P
        imp_3 = 2 + 3 * (nt_1 - 1) - 1  # S
        imp_5 = 2 + 3 * (nt_1 - 1) + 1  # P
    if mp_2 == 'S':
        imp_2 = 2 + 3 * (nt_2 - 1) - 1  # S
        imp_4 = 2 + 3 * (nt_2 - 1) + 1  # P
        imp_6 = 2 + 3 * (nt_2 - 1) + 2  # S
    elif mp_2 == 'B':
        imp_2 = 2 + 3 * (nt_2 - 1)      # B
        imp_4 = 2 + 3 * (nt_2 - 1) - 1  # S
        imp_6 = 2 + 3 * (nt_2 - 1) + 1  # P
    elif mp_2 == 'P':
        imp_2 = 2 + 3 * (nt_2 - 1) - 2  # P
        imp_4 = 2 + 3 * (nt_2 - 1) - 1  # S
        imp_6 = 2 + 3 * (nt_2 - 1) + 1  # P

    if c[0] == 'CAN':  ## Canonical base pairs (A-form RNA)
        (dist_native, ang1_native, ang2_native, 
         dih0_native, dih1_native, dih2_native, nHB) = hb_ARNA_native(nt_1, nt_2)

    elif c[0] == 'NON':  ## Other hydrogen bonds (including non-canonical basepairs)
        xyz1 = chain.get_atom( imp_1 - 1 ).xyz
        xyz2 = chain.get_atom( imp_2 - 1 ).xyz
        xyz3 = chain.get_atom( imp_3 - 1 ).xyz
        xyz4 = chain.get_atom( imp_4 - 1 ).xyz
        xyz5 = chain.get_atom( imp_5 - 1 ).xyz
        xyz6 = chain.get_atom( imp_6 - 1 ).xyz

        v12 = xyz1 - xyz2
        v13 = xyz1 - xyz3
        v53 = xyz5 - xyz3
        v42 = xyz4 - xyz2
        v46 = xyz4 - xyz6
    
        d1212 = v12.dot(v12)
        dist_native = math.sqrt( d1212 )
    
        cos_theta = v13.dot(v12) / math.sqrt(v13.dot(v13) * d1212)
        ang1_native = math.acos(cos_theta) / math.pi * 180.0
    
        cos_theta = v12.dot(v42) / math.sqrt(v42.dot(v42) * d1212)
        ang2_native = math.acos(cos_theta) / math.pi * 180.0
    
        c4212 = v42.cross(v12)
        c1213 = v12.cross(v13)
        dih = math.atan2( v42.dot(c1213) * math.sqrt( v12.dot(v12)), c4212.dot(c1213))
        dih0_native = dih / math.pi * 180.0

        m = v53.cross(v13)
        n = c1213 * -1
        dih = math.atan2( v53.dot(n) * math.sqrt( v13.dot(v13)), m.dot(n))
        dih1_native = dih / math.pi * 180.0

        m = c4212 * -1
        n = v42.cross(v46)
        dih = math.atan2( v12.dot(n) * math.sqrt( v42.dot(v42)), m.dot(n))
        dih2_native = dih / math.pi * 180.0

    hb = HBondDT(iunit1=1,iunit2=1,imp1=imp_1,imp2=imp_2,imp1un=imp_1,imp2un=imp_2,
                 native=dist_native, factor=DT13.HB_U0*nHB, correct_mgo=1.0,coef=DT13.HB_DIST,

                 ang1_imp1=imp_3, ang1_imp2=imp_1, ang1_imp3=imp_2, ang1_iunit1=1,ang1_iunit2=1,
                 ang1_imp1un=imp_3, ang1_imp2un=imp_1, ang1_imp3un=imp_2, 
                 ang1_native=ang1_native, ang1_coef=DT13.HB_ANGL,

                 ang2_imp1=imp_4, ang2_imp2=imp_2, ang2_imp3=imp_1, ang2_iunit1=1,ang2_iunit2=1,
                 ang2_imp1un=imp_4, ang2_imp2un=imp_2, ang2_imp3un=imp_1, 
                 ang2_native=ang2_native, ang2_coef=DT13.HB_ANGL,

                 dih0_imp1=imp_3, dih0_imp2=imp_1, dih0_imp3=imp_2, dih0_imp4=imp_4,dih0_iunit1=1,dih0_iunit2=1,
                 dih0_imp1un=imp_3, dih0_imp2un=imp_1, dih0_imp3un=imp_2, dih0_imp4un=imp_4,
                 dih0_native=dih0_native,dih0_coef=DT13.HB_DIH_HBOND,

                 dih1_imp1=imp_2, dih1_imp2=imp_1, dih1_imp3=imp_3, dih1_imp4=imp_5,dih1_iunit1=1,dih1_iunit2=1,
                 dih1_imp1un=imp_2, dih1_imp2un=imp_1, dih1_imp3un=imp_3, dih1_imp4un=imp_5,
                 dih1_native=dih1_native,dih1_coef=DT13.HB_DIH_CHAIN,

                 dih2_imp1=imp_1, dih2_imp2=imp_2, dih2_imp3=imp_4, dih2_imp4=imp_6,dih2_iunit1=1,dih2_iunit2=1,
                 dih2_imp1un=imp_1, dih2_imp2un=imp_2, dih2_imp3un=imp_4, dih2_imp4un=imp_6,
                 dih2_native=dih2_native,dih2_coef=DT13.HB_DIH_CHAIN,
                 )

    ns.hbondDTs.append(hb)

#nf = NinfoFile('ninfo_RNA13_pdb.ninfo')
nf = NinfoFile(sys.argv[-1])
nf.open_to_write()
nf.write_all(ns)
nf.close()

