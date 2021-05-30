#!/usr/bin/env python

'''
Currently support:
    DT13 from PDB, circ
    NHT19 from PDB

To do:
    NHT19 + circ
    CHT18(DNA)
'''

import sys
import math
import argparse

from cafysis.file_io.pdb import PdbFile
from cafysis.file_io.ninfo import NinfoFile
from cafysis.para.rnaAform import ARNA
from cafysis.para.rnaDT13 import DT13
from cafysis.para.rnaDT15 import DT15
from cafysis.para.rnaNHT19 import NHT19
from cafysis.elements.ninfo import NinfoSet, BondLength, BondAngle, BaseStackDT, HBondDT, TertiaryStackDT


#######################################
## Functions to return native values ##
#######################################

################ Bond
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

################ Angle
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

################ Stack
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

################ H-bond
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
    elif seq[i-1] == 'G' and seq[j-1] == 'U':
        dist = ARNA.HB_GU
        ang1 = ARNA.HBA_SGU
        ang2 = ARNA.HBA_SUG
        dih0 = ARNA.HBD_SGUS
        dih1 = ARNA.HBD_PSGU
        dih2 = ARNA.HBD_PSUG
        nHB = 2
    elif seq[i-1] == 'U' and seq[j-1] == 'G':
        dist = ARNA.HB_GU
        ang1 = ARNA.HBA_SUG
        ang2 = ARNA.HBA_SGU
        dih0 = ARNA.HBD_SGUS
        dih1 = ARNA.HBD_PSUG
        dih2 = ARNA.HBD_PSGU
        nHB = 2
    else:
        print('Canonical basepair should be A-U, G-C, or G-U: ',i,j)
        sys.exit(2)
    return dist, ang1, ang2, dih0, dih1, dih2, nHB


################################################################################
##                                  Main                                      ##
################################################################################
if __name__ == "__main__":

    ###########################################
    ## Define the Parser and check errors    ##
    ###########################################
    parser = argparse.ArgumentParser(
             description='Construct a topology file for TIS model simulations',
             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('outfile', default='ninfo.out', help='output filename')

    group_seq = parser.add_mutually_exclusive_group(required=True)
    group_seq.add_argument('--pdb', type=PdbFile, help='PDB file')
    group_seq.add_argument('--seq', help='Sequence')
    group_seq.add_argument('--seqfile', type=argparse.FileType('r'), help='Sequence FASTA file')

    parser.add_argument('--model', required=True, help='model (DT13, NHT19)')

    parser.add_argument('--hbfile', type=argparse.FileType('r'), help='HB list file')

    parser.add_argument('--tstfile', type=argparse.FileType('r'), help='Tertiary stack list file')

    parser.add_argument('--circ', action="store_true", help='Flag for circRNA')
    
    parser.add_argument('--end5', default='S', help="5'-end P or S")
    parser.add_argument('--end3', default='B', help="3'-end P, S or B")

    parser.add_argument('--nnhb', type=int, default=0, help='Treatment of non-native HB interactions: 0:nothing, 1:NNHB if CAN does not exist, 2:NNHB if CAN or NON do not exist between the bases, 3:NNHB if CAN or NON do not exist between the nucleotides involving one of the bases.')
    parser.add_argument('--nnhb_exclude_file', type=argparse.FileType('r'), help='File that describes NNHB exclusion')

    parser.add_argument('--exvfile', type=argparse.FileType('w'), help='output exv pair file')
    parser.add_argument('--psffile', type=argparse.FileType('w'), help='output PSF file')

    args = parser.parse_args()

    if args.model == 'DT13':
        model = DT13
    elif args.model == 'NHT19':
        model = NHT19
    else:
        print ('Error: model has to be either DT13 or NHT19')
        sys.exit(2)

    if args.pdb is not None:
        #pdb = PdbFile(sys.argv[1])
        args.pdb.open_to_read()
        chains = args.pdb.read_all()
        args.pdb.close()
        
        chain = chains[0]
        n_nt = chain.num_res()
        
        seq = []
        for r in chain.residues:
            # "RA " ---> "A"
            seq.append(r.atoms[0].res_name.strip()[1])

    elif args.seq is not None:
        seq = args.seq
        n_nt = len(seq)

    elif args.seqfile is not None:
        seq = ''
        for l in args.seqfile:
            if l[0] == '>' or l[0] == '#' or l[0] == ';':
                continue
            seq += l.strip()
        n_nt = len(seq)

    else:
        print ('Error: this error should be detected by mutually_exclusive_group of the parser')
        sys.exit(2)

    if args.circ:
        if args.end5 == 'S':
            print ('end5 option is ignored for circRNA')
            args.end5 = 'P'
        if args.end3 in ('P', 'S'):
            print ('end3 option is ignored for circRNA')
            args.end3 = 'B'

    if args.end5 == 'P':
        offset_end5 = 1
    elif args.end5 == 'S':
        offset_end5 = 0
    else:
        print ('end5 option has to be Either S or P')
        sys.exit(2)

    if args.end3 not in ('P', 'S', 'B'):
        print ('end3 option has to be Either P, S, or B')
        sys.exit(2)

    if args.nnhb not in (0,1,2,3):
        print ('Error: --nnhb must be either of ')
        print ('        0: No non-native HB (NNHB).')
        print ('        1: Include NNHBs where no CAN (canonical base pair) interaction exists between the nucleoties.')
        print ('        2: Include NNHBs where neither CAN or NON (non-canonical HB) interaction exist between the bases')
        print ('        3: Include NNHBs where neither CAN or NON (non-canonical HB) interaction exist between the nucleoties involving one of the bases.')
        print ('')
        print ('        The standard selection is 1 for DT15, and 3 for NHT19.')
        print ('')
        sys.exit(2)


    ###########################
    ## Print the setting     ##
    ###########################
    print('#nt: ', n_nt)
    print('Sequence:')
    print(seq)

    if args.circ:
        print ('This is a circRNA. (--end5 option will be ignored)')
    else:
        print("This is a linear RNA with the 5'-end starting with %s" % (args.end5,))

    if args.exvfile:
        print('A file for exv pair list will be generated.')
    
    if args.psffile:
        print('A PSF file will be generated.')

    ################################################
    ## This is the main object to be constructed  ##
    ################################################
    ns = NinfoSet()
    

    ##############
    ## Bond     ##
    ##############
    # for the first nt
    if args.end5 == 'P':
        # P - S
        bl = BondLength(iunit1=1, iunit2=1,
                        imp1=0+offset_end5, imp2=1+offset_end5, 
                        imp1un=0+offset_end5, imp2un=1+offset_end5,
                        native=ARNA.BL_PS, factor=1.0, correct_mgo=1.0, coef=model.BL_PS, type_str='PS')
        ns.bondlengths.append(bl)

    # S - B
    native, type_str = bl_native(1)
    bl = BondLength(iunit1=1, iunit2=1,
                    imp1=1+offset_end5, imp2=2+offset_end5, 
                    imp1un=1+offset_end5, imp2un=2+offset_end5,
                    native=native, factor=1.0, correct_mgo=1.0, coef=model.BL_SB, type_str=type_str)
    ns.bondlengths.append(bl)
    
    # for the second through last nt
    for i in range(2,n_nt+1):
        imp_P = 3 * (i-1) + offset_end5
        imp_S0 = imp_P - 2
        imp_S = imp_P + 1
        imp_B = imp_S + 1
        # S0 - P (from previous nt)
        bl = BondLength(iunit1=1,iunit2=1,imp1=imp_S0,imp2=imp_P,imp1un=imp_S0,imp2un=imp_P,
                        native=ARNA.BL_SP,factor=1.0,correct_mgo=1.0,coef=model.BL_SP,type_str='SP')
        ns.bondlengths.append(bl)
    
        if i == n_nt and args.end3 == 'P':
            break

        # P - S
        bl = BondLength(iunit1=1,iunit2=1,imp1=imp_P,imp2=imp_S,imp1un=imp_P,imp2un=imp_S,
                        native=ARNA.BL_PS,factor=1.0,correct_mgo=1.0,coef=model.BL_PS,type_str='PS')
        ns.bondlengths.append(bl)

        if i == n_nt and args.end3 == 'S':
            break
    
        # S - B
        native, type_str = bl_native(i)
        bl = BondLength(iunit1=1,iunit2=1,imp1=imp_S,imp2=imp_B,imp1un=imp_S,imp2un=imp_B,
                        native=native,factor=1.0,correct_mgo=1.0,coef=model.BL_SB,type_str=type_str)
        ns.bondlengths.append(bl)

    ## circRNA: connect the 5' and 3' ends
    if args.circ:
        imp_S0 = 3*n_nt - 1 # S of 3'-end
        imp_P  = 1          # P between 3'-end and 5'-end

        # S0 - P  (S of 3'-end to P of 5'-end)
        bl = BondLength(iunit1=1,iunit2=1,imp1=imp_S0,imp2=imp_P,imp1un=imp_S0,imp2un=imp_P,
                        native=ARNA.BL_SP,factor=1.0,correct_mgo=1.0,coef=model.BL_SP,type_str='SP')
        ns.bondlengths.append(bl)
    
    
    ##############
    ## Angle    ##
    ##############
    for i in range(1, n_nt+1):
        imp_P  = 3*(i-1) + offset_end5   # This can be 0 if i == 1 and end5 == 'S'
        imp_S  = imp_P + 1
        imp_B  = imp_P + 2
        imp_P1 = imp_P + 3
        imp_S1 = imp_P + 4

        if i == n_nt:
            if args.circ:
                imp_P1 = 1
                imp_S1 = 2

            else:
                """ 
                If the last nucleotides does not have B, then angles formed by P and S
                are already taken care in the previous (n_nt - 1), thus break the loop.
                """
                if args.end3 != 'B':
                    break

        """ P - S - B
        1st nucleotide: When there is P
        2nd through 2nd last: Always 
        Last: Only when there is B (including circular)
        """
        if (1 < i < n_nt) or (i == 1 and args.end5 == 'P') or (i == n_nt and args.end3 == 'B'):
            native, type_str = ba_PSB_native(i)
            ba = BondAngle(iunit1=1, iunit2=1, imp1=imp_P, imp2=imp_S, imp3=imp_B, 
                           imp1un=imp_P, imp2un=imp_S, imp3un=imp_B, 
                           native=native, factor=1.0, correct_mgo=1.0, 
                           coef=model.BA_PSB, type_str=type_str)
            ns.bondangles.append(ba)

        """ P - S - P1
        1st nucleotide: When there is P
        2nd through 2nd last: Always 
        Last nucleotide: Only when circular
        """
        if (1 < i < n_nt) or (i == 1 and args.end5 == 'P') or (i == n_nt and args.circ):
            ba = BondAngle(iunit1=1, iunit2=1, imp1=imp_P, imp2=imp_S, imp3=imp_P1, 
                           imp1un=imp_P, imp2un=imp_S, imp3un=imp_P1, 
                           native=ARNA.BA_PSP, factor=1.0, correct_mgo=1.0, 
                           coef=model.BA_PSP, type_str='PSP')
            ns.bondangles.append(ba)

        """ B - S - P1
        1st through 2nd last: Always 
        Last nucleotide: Only when circular
        """
        if (i < n_nt) or (i == n_nt and args.circ):
            native, type_str = ba_BSP_native(i)
            ba = BondAngle(iunit1=1, iunit2=1, imp1=imp_B, imp2=imp_S, imp3=imp_P1,
                           imp1un=imp_B, imp2un=imp_S, imp3un=imp_P1,
                           native=native,factor=1.0, correct_mgo=1.0,
                           coef=model.BA_BSP, type_str=type_str)
            ns.bondangles.append(ba)

        """ S - P1 - S1
        1st through 3rd last: Always
        2nd last: When there is S1 (including circular)
        Last nucleotide: Only when circular
        """
        if (i < n_nt - 1) or (i == n_nt-1 and args.end3 != 'P') or (i == n_nt and args.circ):
            ba = BondAngle(iunit1=1, iunit2=1, imp1=imp_S, imp2=imp_P1, imp3=imp_S1,
                           imp1un=imp_S, imp2un=imp_P1, imp3un=imp_S1,
                           native=ARNA.BA_SPS, factor=1.0, correct_mgo=1.0,
                           coef=model.BA_SPS, type_str='SPS')
            ns.bondangles.append(ba)

    
    ##############
    ## Stack    ##
    ##############
    # for the second through one before last nt
    nt_start = 3  # end5 = 'S'
    if args.end5 == 'P':
        nt_start = 2

    for i in range(nt_start,n_nt):
        imp_P2 = 3*(i-1) + offset_end5
        imp_P1 = imp_P2 - 3
        imp_S1 = imp_P2 - 2
        imp_B1 = imp_P2 - 1
        imp_S2 = imp_P2 + 1
        imp_B2 = imp_P2 + 2
        imp_P3 = imp_P2 + 3
    
        native, type_str = bs_native(i-1,i)
        h, s, Tm = model.ST_U0[ type_str[0]+type_str[2] ]
        bs = BaseStackDT(iunit1=1,iunit2=1,imp1=imp_B1,imp2=imp_B2,imp1un=imp_B1,imp2un=imp_B2,
                           native=native, factor=0.0,correct_mgo=1.0,coef=model.ST_DIST,type_str=type_str,
                           h=h, s=s, Tm=Tm,
                           dih1_imp1=imp_P1, dih1_imp2=imp_S1, dih1_imp3=imp_P2, dih1_imp4=imp_S2,dih1_iunit1=1,dih1_iunit2=1,
                           dih1_imp1un=imp_P1, dih1_imp2un=imp_S1, dih1_imp3un=imp_P2, dih1_imp4un=imp_S2,
                           dih1_native=ARNA.DIH_PSPS,dih1_coef=model.ST_DIH,dih1_type_str='PSPS',
                           dih2_imp1=imp_S1, dih2_imp2=imp_P2, dih2_imp3=imp_S2, dih2_imp4=imp_P3,dih2_iunit1=1,dih2_iunit2=1,
                           dih2_imp1un=imp_S1, dih2_imp2un=imp_P2, dih2_imp3un=imp_S2, dih2_imp4un=imp_P3,
                           dih2_native=ARNA.DIH_SPSP,dih2_coef=model.ST_DIH,dih2_type_str='SPSP')
        ns.basestackDTs.append(bs)
    
    # circRNA
    if args.circ:
        # (n_nt - 1) : (n_nt)
        imp_P1 = 3*n_nt - 5
        imp_S1 = 3*n_nt - 4
        imp_B1 = 3*n_nt - 3
        imp_P2 = 3*n_nt - 2
        imp_S2 = 3*n_nt - 1
        imp_B2 = 3*n_nt
        imp_P3 = 1
        native, type_str = bs_native(n_nt-1,n_nt)
        h, s, Tm = model.ST_U0[ type_str[0]+type_str[2] ]
        bs = BaseStackDT(iunit1=1,iunit2=1,imp1=imp_B1,imp2=imp_B2,imp1un=imp_B1,imp2un=imp_B2,
                           native=native, factor=0.0,correct_mgo=1.0,coef=model.ST_DIST,type_str=type_str,
                           h=h, s=s, Tm=Tm,
                           dih1_imp1=imp_P1, dih1_imp2=imp_S1, dih1_imp3=imp_P2, dih1_imp4=imp_S2,dih1_iunit1=1,dih1_iunit2=1,
                           dih1_imp1un=imp_P1, dih1_imp2un=imp_S1, dih1_imp3un=imp_P2, dih1_imp4un=imp_S2,
                           dih1_native=ARNA.DIH_PSPS,dih1_coef=model.ST_DIH,dih1_type_str='PSPS',
                           dih2_imp1=imp_S1, dih2_imp2=imp_P2, dih2_imp3=imp_S2, dih2_imp4=imp_P3,dih2_iunit1=1,dih2_iunit2=1,
                           dih2_imp1un=imp_S1, dih2_imp2un=imp_P2, dih2_imp3un=imp_S2, dih2_imp4un=imp_P3,
                           dih2_native=ARNA.DIH_SPSP,dih2_coef=model.ST_DIH,dih2_type_str='SPSP')
        ns.basestackDTs.append(bs)

        # (n_nt) : (1)
        imp_P1 = 3*n_nt - 2
        imp_S1 = 3*n_nt - 1
        imp_B1 = 3*n_nt
        imp_P2 = 1
        imp_S2 = 2
        imp_B2 = 3
        imp_P3 = 4
        native, type_str = bs_native(n_nt,1)
        h, s, Tm = model.ST_U0[ type_str[0]+type_str[2] ]
        bs = BaseStackDT(iunit1=1,iunit2=1,imp1=imp_B1,imp2=imp_B2,imp1un=imp_B1,imp2un=imp_B2,
                           native=native, factor=0.0,correct_mgo=1.0,coef=model.ST_DIST,type_str=type_str,
                           h=h, s=s, Tm=Tm,
                           dih1_imp1=imp_P1, dih1_imp2=imp_S1, dih1_imp3=imp_P2, dih1_imp4=imp_S2,dih1_iunit1=1,dih1_iunit2=1,
                           dih1_imp1un=imp_P1, dih1_imp2un=imp_S1, dih1_imp3un=imp_P2, dih1_imp4un=imp_S2,
                           dih1_native=ARNA.DIH_PSPS,dih1_coef=model.ST_DIH,dih1_type_str='PSPS',
                           dih2_imp1=imp_S1, dih2_imp2=imp_P2, dih2_imp3=imp_S2, dih2_imp4=imp_P3,dih2_iunit1=1,dih2_iunit2=1,
                           dih2_imp1un=imp_S1, dih2_imp2un=imp_P2, dih2_imp3un=imp_S2, dih2_imp4un=imp_P3,
                           dih2_native=ARNA.DIH_SPSP,dih2_coef=model.ST_DIH,dih2_type_str='SPSP')
        ns.basestackDTs.append(bs)


    
    ##############
    ## H-bond   ##
    ##############
    hblist = []
    if args.hbfile is not None:
        for l in args.hbfile:
            if l.find('#') != -1:
                continue
            lsp = l.strip().split()
        
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

            atoms1 = []
            atoms2 = []

            # Get pairs of atom names, two at a time, and append to atoms1 and atoms2.
            it = iter(lsp[6:])
            for a1, a2 in zip(it,it):
                atoms1.append(a1)
                atoms2.append(a2)

            '''
            !! atoms1 and atoms2 can be empty if model is DT13. !!
            '''
        
            # Order the nucleotide numbers
            nt1 = int(lsp[1])
            nt2 = int(lsp[3])
            if nt1 <= nt2:
                hblist.append((lsp[0], nt1, lsp[2], nt2, lsp[4], int(lsp[5]), atoms1, atoms2))
            else:
                hblist.append((lsp[0], nt2, lsp[4], nt1, lsp[2], int(lsp[5]), atoms2, atoms1))
        

    #### Non-native HB
    nnlist = []
    if args.nnhb in (1,2,3):

        nnhb_exclude_nts = set()
        nnhb_exclude_pairs = set()
        if args.nnhb_exclude_file is not None:

            for l in args.nnhb_exclude_file:

                if l.find('#') != -1:
                    continue

                lsp = l.strip().split()

                if len(lsp) == 2:
                    c1  = int(lsp[0])
                    nt1 = int(lsp[1])

                    nnhb_exclude_nts.add((c1, nt1))

                elif len(lsp) == 4:
                    c1  = int(lsp[0])
                    nt1 = int(lsp[1])
                    c2  = int(lsp[2])
                    nt2 = int(lsp[3])
        
                    nnhb_exclude_nts.add((c1, nt1, c2, nt2))

                else:
                    print('Error: unknown format in nnhb_exclude_file:\n' + l)
                    sys.exit(2)

        if args.circ:
            print("Error: The combination of NN and circRNA is not implemented yet.")
            sys.exit(2)

        else:
            '''
            HB interactions between bases need Ps in the downstream.
            '''
            nt1_begin = 1
            nt1_end = n_nt - model.NNHB_NT_SEP - 1
            nt2_end = n_nt - 1

        for nt1 in range(nt1_begin, nt1_end+1):

            if (1,nt1) in nnhb_exclude_nts:
                continue
            
            for nt2 in range(nt1 + model.NNHB_NT_SEP, nt2_end + 1):
    
                if (1,nt1,1,nt2) in nnhb_exclude_pairs:
                    continue

                if   seq[nt1-1] == 'G' and seq[nt2-1] == 'C':
                    atoms1 = ("N1", "N2", "O6")
                    atoms2 = ("N3", "O2", "N4")
                    nHB = 3
                elif seq[nt1-1] == 'C' and seq[nt2-1] == 'G':
                    atoms1 = ("N3", "O2", "N4")
                    atoms2 = ("N1", "N2", "O6")
                    nHB = 3
                elif seq[nt1-1] == 'A' and seq[nt2-1] == 'U':
                    atoms1 = ("N1", "N6")
                    atoms2 = ("N3", "O4")
                    nHB = 2
                elif seq[nt1-1] == 'U' and seq[nt2-1] == 'A':
                    atoms1 = ("N3", "O4")
                    atoms2 = ("N1", "N6")
                    nHB = 2
                elif seq[nt1-1] == 'G' and seq[nt2-1] == 'U':
                    atoms1 = ("N1", "O6")
                    atoms2 = ("O2", "N3")
                    nHB = 2
                elif seq[nt1-1] == 'U' and seq[nt2-1] == 'G':
                    atoms1 = ("O2", "N3")
                    atoms2 = ("N1", "O6")
                    nHB = 2
                else:
                    continue
    
                flg_exist = False
                for c in hblist:
                    if args.nnhb == 1:
                        if c[0] == 'CAN' and c[1] == nt1 and c[3] == nt2 and c[2] == 'B' and c[4] == 'B':
                            flg_exist = True
                            break

                    elif args.nnhb == 2:
                        if (c[0] == 'CAN' or c[0] == 'NON') and c[1] == nt1 and c[3] == nt2 and c[2] == 'B' and c[4] == 'B':
                            flg_exist = True
                            break

                    elif args.nnhb == 3:
                        if (c[0] == 'CAN' or c[0] == 'NON') and c[1] == nt1 and c[3] == nt2 and c[2] == 'B' and c[4] == 'B':
                            flg_exist = True
                            break
                        if (c[0] == 'NON') and c[1] == nt1 and c[3] == nt2 and (c[2] == 'B' or c[4] == 'B'):
                            flg_exist = True
                            break

                    else:
                        print ("Error: logical error in hblist loop in nnhb. Examine the code.")
                        sys.exit(2)
    
                if flg_exist:
                    continue
    
                nnlist.append(('CAN', nt1, 'B', nt2, 'B', nHB, atoms1, atoms2))

    for c in hblist + nnlist:
        nt_1 = c[1]
        mp_1 = c[2]
        nt_2 = c[3]
        mp_2 = c[4]
        atoms1 = c[6]
        atoms2 = c[7]
    
        # For both CAN and NON, 
        if mp_1 == 'S':
            imp_1 = 2 + 3 * (nt_1 - 1) - 1 + offset_end5 # S
            imp_3 = 2 + 3 * (nt_1 - 1) + 1 + offset_end5 # P
            imp_5 = 2 + 3 * (nt_1 - 1) + 2 + offset_end5 # S
        elif mp_1 == 'B':
            imp_1 = 2 + 3 * (nt_1 - 1)     + offset_end5 # B
            imp_3 = 2 + 3 * (nt_1 - 1) - 1 + offset_end5 # S
            imp_5 = 2 + 3 * (nt_1 - 1) + 1 + offset_end5 # P
        elif mp_1 == 'P':
            imp_1 = 2 + 3 * (nt_1 - 1) - 2 + offset_end5 # P
            imp_3 = 2 + 3 * (nt_1 - 1) - 1 + offset_end5 # S
            imp_5 = 2 + 3 * (nt_1 - 1) + 1 + offset_end5 # P
        if mp_2 == 'S':
            imp_2 = 2 + 3 * (nt_2 - 1) - 1 + offset_end5 # S
            imp_4 = 2 + 3 * (nt_2 - 1) + 1 + offset_end5 # P
            imp_6 = 2 + 3 * (nt_2 - 1) + 2 + offset_end5 # S
        elif mp_2 == 'B':
            imp_2 = 2 + 3 * (nt_2 - 1)     + offset_end5 # B
            imp_4 = 2 + 3 * (nt_2 - 1) - 1 + offset_end5 # S
            imp_6 = 2 + 3 * (nt_2 - 1) + 1 + offset_end5 # P
        elif mp_2 == 'P':
            imp_2 = 2 + 3 * (nt_2 - 1) - 2 + offset_end5 # P
            imp_4 = 2 + 3 * (nt_2 - 1) - 1 + offset_end5 # S
            imp_6 = 2 + 3 * (nt_2 - 1) + 1 + offset_end5 # P
    
        if c[0] == 'CAN':  ## Canonical base pairs (A-form RNA)
            (dist_native, ang1_native, ang2_native, 
             dih0_native, dih1_native, dih2_native, nHB) = hb_ARNA_native(nt_1, nt_2)

            if nHB != c[5]:
                print('Error: nHB != c[5]')
                print(c)
                sys.exit(2)

            sectert = 'S'
    
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

            sectert = 'T'
            nHB  = c[5]

        factor = model.HB_U0 * nHB

        if args.model == 'DT13':
            # Fot DT13, do not output sectert, nHB, atoms1, and atoms2.
            nHB = None
            sectert = None
            atoms1 = []
            atoms2 = []

        hb = HBondDT(iunit1=1,iunit2=1,imp1=imp_1,imp2=imp_2,imp1un=imp_1,imp2un=imp_2,
                     native=dist_native, factor=factor, correct_mgo=1.0, coef=model.HB_DIST,
    
                     ang1_imp1=imp_3, ang1_imp2=imp_1, ang1_imp3=imp_2, ang1_iunit1=1,ang1_iunit2=1,
                     ang1_imp1un=imp_3, ang1_imp2un=imp_1, ang1_imp3un=imp_2, 
                     ang1_native=ang1_native, ang1_coef=model.HB_ANGL,
    
                     ang2_imp1=imp_4, ang2_imp2=imp_2, ang2_imp3=imp_1, ang2_iunit1=1,ang2_iunit2=1,
                     ang2_imp1un=imp_4, ang2_imp2un=imp_2, ang2_imp3un=imp_1, 
                     ang2_native=ang2_native, ang2_coef=model.HB_ANGL,
    
                     dih0_imp1=imp_3, dih0_imp2=imp_1, dih0_imp3=imp_2, dih0_imp4=imp_4,dih0_iunit1=1,dih0_iunit2=1,
                     dih0_imp1un=imp_3, dih0_imp2un=imp_1, dih0_imp3un=imp_2, dih0_imp4un=imp_4,
                     dih0_native=dih0_native,dih0_coef=model.HB_DIH_HBOND,
    
                     dih1_imp1=imp_2, dih1_imp2=imp_1, dih1_imp3=imp_3, dih1_imp4=imp_5,dih1_iunit1=1,dih1_iunit2=1,
                     dih1_imp1un=imp_2, dih1_imp2un=imp_1, dih1_imp3un=imp_3, dih1_imp4un=imp_5,
                     dih1_native=dih1_native,dih1_coef=model.HB_DIH_CHAIN,
    
                     dih2_imp1=imp_1, dih2_imp2=imp_2, dih2_imp3=imp_4, dih2_imp4=imp_6,dih2_iunit1=1,dih2_iunit2=1,
                     dih2_imp1un=imp_1, dih2_imp2un=imp_2, dih2_imp3un=imp_4, dih2_imp4un=imp_6,
                     dih2_native=dih2_native,dih2_coef=model.HB_DIH_CHAIN,

                     sectert=sectert, nHB=nHB, atoms1=atoms1, atoms2=atoms2
                     )
    
        ns.hbondDTs.append(hb)


    ######################
    ## Tertiary stack   ##
    ######################
    if args.tstfile is not None:

        sstlist = []  # Secondary
        tstlist = []  # Tertiary
        tstlist_force = []  # Tertiary (force to be included)

        for l in args.tstfile:

            if l.find('#') != -1:
                continue
            lsp = l.split()
        
            nt1 = int(lsp[1])
            nt2 = int(lsp[2])
        
            if abs(nt2) < abs(nt1):
                save = nt1
                nt1 = nt2
                nt2 = save
        
            if lsp[0] == 'S':
                sstlist.append( (nt1, nt2) )
            elif lsp[0] == 'T':
                tstlist.append( (nt1, nt2) )
            elif lsp[0] == 'F':
                tstlist_force.append( (nt1, nt2) )
            else:
                print('Error: unknown stack type')
                sys.exit(2)
        
        #### Remove inner stem
        del_id = []
        for itst, tst in enumerate(tstlist):
            nt1 = tst[0]
            nt2 = tst[1]
        
            if nt1 < 0:
                sgn1 = -1
                nt1 = -nt1
            else:
                sgn1 = 1
            if nt2 < 0:
                sgn2 = -1
                nt2 = -nt2
            else:
                sgn2 = 1
        
            mp1 = 2 + 3 * (nt1 - 1) + offset_end5      # B
            mp2 = 2 + 3 * (nt2 - 1) + offset_end5      # B
        
            if sgn1 < 0:
                nei_mp1 = mp1 - 3
            else:
                nei_mp1 = mp1 + 3
            if sgn2 < 0:
                nei_mp2 = mp2 - 3
            else:
                nei_mp2 = mp2 + 3
        
            if nei_mp1 <= 0 or nei_mp2 <= 0:
                continue
        
            flg_0 = False
            flg_1 = False
            flg_2 = False
            for hb in ns.hbondDTs:
                if hb.sectert != 'S':
                    continue
                if hb.imp1 == mp1 and hb.imp2 == mp2:
                    flg_0 = True
                elif hb.imp1 == mp1 and hb.imp2 == nei_mp2:
                    flg_1 = True
                elif hb.imp1 == nei_mp1 and hb.imp2 == mp2:
                    flg_2 = True
                if flg_0 or (flg_1 and flg_2):
                    break
        
            if flg_0 or (flg_1 and flg_2):
                del_id.append(itst)
        
        for itst in sorted(del_id, reverse=True):
            del tstlist[itst]
        
        excess = []
        for i in range(n_nt+1):
            excess.append([0,0])   

        for isst, sst in enumerate(sstlist):
            nt1 = sst[0]
            nt2 = sst[1]
        
            if nt1 < 0:
                sgn1 = -1
                nt1 = -nt1
                excess[nt1][0] += 1
            else:
                sgn1 = 1
                excess[nt1][1] += 1
            if nt2 < 0:
                sgn2 = -1
                nt2 = -nt2
                excess[nt2][0] += 1
            else:
                sgn2 = 1
                excess[nt2][1] += 1

        for itst, tst in enumerate(tstlist):
            nt1 = tst[0]
            nt2 = tst[1]
        
            if nt1 < 0:
                sgn1 = -1
                nt1 = -nt1
                excess[nt1][0] += 1
            else:
                sgn1 = 1
                excess[nt1][1] += 1
            if nt2 < 0:
                sgn2 = -1
                nt2 = -nt2
                excess[nt2][0] += 1
            else:
                sgn2 = 1
                excess[nt2][1] += 1

        for itst, tst in enumerate(tstlist_force):
            nt1 = tst[0]
            nt2 = tst[1]
        
            if nt1 < 0:
                sgn1 = -1
                nt1 = -nt1
                excess[nt1][0] += 1
            else:
                sgn1 = 1
                excess[nt1][1] += 1
            if nt2 < 0:
                sgn2 = -1
                nt2 = -nt2
                excess[nt2][0] += 1
            else:
                sgn2 = 1
                excess[nt2][1] += 1
        
        # If excess >  1 ---> inclusive
        #    excess <= 1 ---> exclusive
        
        for itst, tst in enumerate(tstlist + tstlist_force):
            nt1 = tst[0]
            nt2 = tst[1]
        
            if nt1 < 0:
                sgn1 = -1
                nt1 = -nt1
            else:
                sgn1 = 1
            if nt2 < 0:
                sgn2 = -1
                nt2 = -nt2
            else:
                sgn2 = 1
        
            imp_1 = 2 + 3 * (nt1 - 1)      + offset_end5   # B
            imp_3 = 2 + 3 * (nt1 - 1) - 1  + offset_end5   # S
            imp_5 = 2 + 3 * (nt1 - 1) + 1  + offset_end5   # P
            imp_2 = 2 + 3 * (nt2 - 1)      + offset_end5   # B
            imp_4 = 2 + 3 * (nt2 - 1) - 1  + offset_end5   # S
            imp_6 = 2 + 3 * (nt2 - 1) + 1  + offset_end5   # P
        
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
        
            if sgn1 < 0:
                if excess[nt1][0] > 1:
                    excess1 = -1
                else:
                    excess1 = -2
            else:
                if excess[nt1][1] > 1:
                    excess1 = +1
                else:
                    excess1 = +2
        
            if sgn2 < 0:
                if excess[nt2][0] > 1:
                    excess2 = -1
                else:
                    excess2 = -2
            else:
                if excess[nt2][1] > 1:
                    excess2 = +1
                else:
                    excess2 = +2
        
            tst = TertiaryStackDT(iunit1=1,iunit2=1,imp1=imp_1,imp2=imp_2,imp1un=imp_1,imp2un=imp_2,
                         native=dist_native, factor=model.TST_U0, correct_mgo=1.0,coef=model.TST_DIST,
        
                         ang1_imp1=imp_3, ang1_imp2=imp_1, ang1_imp3=imp_2, ang1_iunit1=1,ang1_iunit2=1,
                         ang1_imp1un=imp_3, ang1_imp2un=imp_1, ang1_imp3un=imp_2, 
                         ang1_native=ang1_native, ang1_coef=model.TST_ANGL,
        
                         ang2_imp1=imp_4, ang2_imp2=imp_2, ang2_imp3=imp_1, ang2_iunit1=1,ang2_iunit2=1,
                         ang2_imp1un=imp_4, ang2_imp2un=imp_2, ang2_imp3un=imp_1, 
                         ang2_native=ang2_native, ang2_coef=model.TST_ANGL,
        
                         dih0_imp1=imp_3, dih0_imp2=imp_1, dih0_imp3=imp_2, dih0_imp4=imp_4,dih0_iunit1=1,dih0_iunit2=1,
                         dih0_imp1un=imp_3, dih0_imp2un=imp_1, dih0_imp3un=imp_2, dih0_imp4un=imp_4,
                         dih0_native=dih0_native,dih0_coef=model.TST_DIH,
        
                         dih1_imp1=imp_2, dih1_imp2=imp_1, dih1_imp3=imp_3, dih1_imp4=imp_5,dih1_iunit1=1,dih1_iunit2=1,
                         dih1_imp1un=imp_2, dih1_imp2un=imp_1, dih1_imp3un=imp_3, dih1_imp4un=imp_5,
                         dih1_native=dih1_native,dih1_coef=model.TST_DIH,
        
                         dih2_imp1=imp_1, dih2_imp2=imp_2, dih2_imp3=imp_4, dih2_imp4=imp_6,dih2_iunit1=1,dih2_iunit2=1,
                         dih2_imp1un=imp_1, dih2_imp2un=imp_2, dih2_imp3un=imp_4, dih2_imp4un=imp_6,
                         dih2_native=dih2_native,dih2_coef=model.TST_DIH,
        
                         excess1=excess1, excess2=excess2, 
                         )
        
            ns.tertiarystackDTs.append(tst)
    

    ###########################
    ## Sort                  ##
    ###########################
    ns.update_info()
    ns.sort_by_mp()
    ns.reassign_id()

    ###########################
    ## Generate ninfo file   ##
    ###########################
    nf = NinfoFile(args.outfile)
    nf.open_to_write()
    nf.write_all(ns)
    nf.close()
    

    ###################################
    ## Generate exv file (optional)  ##
    ###################################
    if args.exvfile is not None:

        if args.end5 == 'P':
            nmp = 3 * n_nt
        else:
            nmp = 3 * n_nt - 1

        imp_plus_sep = {}  # An array to store (imp + n_sep_nlocal_X).
                           # X depends on the type of imp.
        for imp in range(1, nmp+1):
            itype = imp % 3

            if args.end5 == 'P':
                if itype == 0: # B
                    imp_plus_sep[imp] = imp + model.n_sep_nlocal_B
                elif itype == 1: # P
                    imp_plus_sep[imp] = imp + model.n_sep_nlocal_P
                else: # S
                    imp_plus_sep[imp] = imp + model.n_sep_nlocal_S
            else:
                if itype == 0: # P
                    imp_plus_sep[imp] = imp + model.n_sep_nlocal_P
                elif itype == 1: # S
                    imp_plus_sep[imp] = imp + model.n_sep_nlocal_S
                else: # B
                    imp_plus_sep[imp] = imp + model.n_sep_nlocal_B

        nexv = 0
        for imp in range(1, nmp):

            for jmp in range(imp, nmp+1):

                if args.circ:
                    if imp_plus_sep[imp] <= jmp and imp_plus_sep[jmp] <= imp + nmp:
                        args.exvfile.write('%i %i\n' % (imp, jmp))
                        nexv += 1
                else:
                    if imp_plus_sep[imp] <= jmp:
                        args.exvfile.write('%i %i\n' % (imp, jmp))
                        nexv += 1

        print ('#exv: %i' % (nexv,))
        args.exvfile.close()


    ###################################
    ## Generate PSF file (optional)  ##
    ###################################
    if args.psffile is not None:

        if args.end5 == 'P':
            nmp = 3 * n_nt
        else:
            nmp = 3 * n_nt - 1

        imp_plus_sep = {}  # An array to store (imp + n_sep_nlocal_X).
                           # X depends on the type of imp.
        for imp in range(1, nmp+1):
            itype = imp % 3

            if args.end5 == 'P':
                if itype == 0: # B
                    imp_plus_sep[imp] = imp + model.n_sep_nlocal_B
                elif itype == 1: # P
                    imp_plus_sep[imp] = imp + model.n_sep_nlocal_P
                else: # S
                    imp_plus_sep[imp] = imp + model.n_sep_nlocal_S
            else:
                if itype == 0: # P
                    imp_plus_sep[imp] = imp + model.n_sep_nlocal_P
                elif itype == 1: # S
                    imp_plus_sep[imp] = imp + model.n_sep_nlocal_S
                else: # B
                    imp_plus_sep[imp] = imp + model.n_sep_nlocal_B

        nexv = 0
        for imp in range(1, nmp):

            for jmp in range(imp, nmp+1):

                if args.circ:
                    if imp_plus_sep[imp] <= jmp and imp_plus_sep[jmp] <= imp + nmp:
                        args.exvfile.write('%i %i\n' % (imp, jmp))
                        nexv += 1
                else:
                    if imp_plus_sep[imp] <= jmp:
                        args.exvfile.write('%i %i\n' % (imp, jmp))
                        nexv += 1

        print ('#exv: %i' % (nexv,))
        args.exvfile.close()
