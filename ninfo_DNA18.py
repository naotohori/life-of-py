#!/usr/bin/env python

import sys
import math
import argparse

from cafysis.file_io.pdb import PdbFile
from cafysis.file_io.ninfo import NinfoFile
from cafysis.para.dnaBform import BDNA
from cafysis.para.dnaCHT18 import CHT18
from cafysis.elements.ninfo import NinfoSet, BondLength, BondAngle, BaseStackDT, HBondDT

NNHB_NT_SEP = 5

#######################################
## Functions to return native values ##
#######################################

################ Bond
def bl_native(iseq, i):

    if seq[iseq][i-1] == 'A':
        native = BDNA.BL_SA
        coef = CHT18.BL_SA
        type_str = 'SA'
    elif seq[iseq][i-1] == 'T':
        native = BDNA.BL_ST
        coef = CHT18.BL_ST
        type_str = 'ST'
    elif seq[iseq][i-1] == 'G':
        native = BDNA.BL_SG
        coef = CHT18.BL_SG
        type_str = 'SG'
    elif seq[iseq][i-1] == 'C':
        native = BDNA.BL_SC
        coef = CHT18.BL_SC
        type_str = 'SC'

    return native, coef, type_str

################ Angle
def ba_BSP_native(iseq, i):

    if seq[iseq][i-1] == 'A':
        native = BDNA.BA_ASP
        coef = CHT18.BA_ASP
        type_str = 'ASP'
    elif seq[iseq][i-1] == 'T':
        native = BDNA.BA_TSP
        coef = CHT18.BA_TSP
        type_str = 'TSP'
    elif seq[iseq][i-1] == 'G':
        native = BDNA.BA_GSP
        coef = CHT18.BA_GSP
        type_str = 'GSP'
    elif seq[iseq][i-1] == 'C':
        native = BDNA.BA_CSP
        coef = CHT18.BA_CSP
        type_str = 'CSP'

    return native, coef, type_str

def ba_PSB_native(iseq, i):

    if seq[iseq][i-1] == 'A':
        native = BDNA.BA_PSA
        coef = CHT18.BA_PSA
        type_str = 'PSA'
    elif seq[iseq][i-1] == 'T':
        native = BDNA.BA_PST
        coef = CHT18.BA_PST
        type_str = 'PST'
    elif seq[iseq][i-1] == 'G':
        native = BDNA.BA_PSG
        coef = CHT18.BA_PSG
        type_str = 'PSG'
    elif seq[iseq][i-1] == 'C':
        native = BDNA.BA_PSC
        coef = CHT18.BA_PSC
        type_str = 'PSC'

    return native, coef, type_str

################ Stack
def bs_native(iseq, i,j):
    if seq[iseq][i-1] == 'A':
        if seq[iseq][j-1] == 'A':
            native = BDNA.ST_AA
            type_str = 'A-A'
        elif seq[iseq][j-1] == 'T':
            native = BDNA.ST_AT
            type_str = 'A-T'
        elif seq[iseq][j-1] == 'G':
            native = BDNA.ST_AG
            type_str = 'A-G'
        elif seq[iseq][j-1] == 'C':
            native = BDNA.ST_AC
            type_str = 'A-C'
    elif seq[iseq][i-1] == 'T':
        if seq[iseq][j-1] == 'A':
            native = BDNA.ST_TA
            type_str = 'T-A'
        elif seq[iseq][j-1] == 'T':
            native = BDNA.ST_TT
            type_str = 'T-T'
        elif seq[iseq][j-1] == 'G':
            native = BDNA.ST_TG
            type_str = 'T-G'
        elif seq[iseq][j-1] == 'C':
            native = BDNA.ST_TC
            type_str = 'T-C'
    elif seq[iseq][i-1] == 'G':
        if seq[iseq][j-1] == 'A':
            native = BDNA.ST_GA
            type_str = 'G-A'
        elif seq[iseq][j-1] == 'T':
            native = BDNA.ST_GT
            type_str = 'G-T'
        elif seq[iseq][j-1] == 'G':
            native = BDNA.ST_GG
            type_str = 'G-G'
        elif seq[iseq][j-1] == 'C':
            native = BDNA.ST_GC
            type_str = 'G-C'
    elif seq[iseq][i-1] == 'C':
        if seq[iseq][j-1] == 'A':
            native = BDNA.ST_CA
            type_str = 'C-A'
        elif seq[iseq][j-1] == 'T':
            native = BDNA.ST_CT
            type_str = 'C-T'
        elif seq[iseq][j-1] == 'G':
            native = BDNA.ST_CG
            type_str = 'C-G'
        elif seq[iseq][j-1] == 'C':
            native = BDNA.ST_CC
            type_str = 'C-C'
    return native, type_str

################ H-bond
def hb_BDNA_native(iseq, i, jseq, j):
    if seq[iseq][i-1] == 'A' and seq[jseq][j-1] == 'T':
        dist = BDNA.HB_AT
        ang1 = BDNA.HBA_SAT
        ang2 = BDNA.HBA_STA
        dih0 = BDNA.HBD_SATS
        dih1 = BDNA.HBD_PSAT
        dih2 = BDNA.HBD_PSTA
        nHB = 2
    elif seq[iseq][i-1] == 'T' and seq[jseq][j-1] == 'A':
        dist = BDNA.HB_AT
        ang1 = BDNA.HBA_STA
        ang2 = BDNA.HBA_SAT
        dih0 = BDNA.HBD_SATS
        dih1 = BDNA.HBD_PSTA
        dih2 = BDNA.HBD_PSAT
        nHB = 2
    elif seq[iseq][i-1] == 'G' and seq[jseq][j-1] == 'C':
        dist = BDNA.HB_GC
        ang1 = BDNA.HBA_SGC
        ang2 = BDNA.HBA_SCG
        dih0 = BDNA.HBD_SGCS
        dih1 = BDNA.HBD_PSGC
        dih2 = BDNA.HBD_PSCG
        nHB = 3
    elif seq[iseq][i-1] == 'C' and seq[jseq][j-1] == 'G':
        dist = BDNA.HB_GC
        ang1 = BDNA.HBA_SCG
        ang2 = BDNA.HBA_SGC
        dih0 = BDNA.HBD_SGCS
        dih1 = BDNA.HBD_PSCG
        dih2 = BDNA.HBD_PSGC
        nHB = 3
    else:
        print('Canonical basepair should be A-T or G-C: ',i,j)
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
#    group_seq.add_argument('--pdb', type=PdbFile, help='PDB file')  # --pdb option not available for DNA now
    group_seq.add_argument('--seq', help='Sequence (single stranded DNA)')
    group_seq.add_argument('--seq_ds', help='Sequence (for dsDNA; automatically generate the complementary sequence)')
    group_seq.add_argument('--seqfile', type=argparse.FileType('r'), help='Sequence FASTA file')

    parser.add_argument('--hbfile', type=argparse.FileType('r'), help='HB list file')
    parser.add_argument('--hb_canonical', default=False, action="store_true", help='All canonical base pairs')

#    parser.add_argument('--circ', action="store_true", help='Flag for circRNA')
    
    parser.add_argument('--end5', default='S', help="5'-end P or S")

    parser.add_argument('--exvfile', type=argparse.FileType('w'), help='output exv pair file')

    args = parser.parse_args()

# Currently circular DNA is not supported
    args.circ = False

#    if args.pdb is not None:
#        #pdb = PdbFile(sys.argv[1])
#        args.pdb.open_to_read()
#        chains = args.pdb.read_all()
#        args.pdb.close()
#        
#        chain = chains[0]
#        n_nt = chain.num_res()
#        
#        seq = []
#        for r in chain.residues:
#            # "RA " ---> "A"
#            seq.append(r.atoms[0].res_name.strip()[1])

#    elif args.seq is not None:
    if args.seq is not None:
        # --seq is only for single stranded DNA simulaitons
        seq = []
        seq.append(args.seq)
        n_nt_all = len(args.seq)
        n_seq = 1

    elif args.seq_ds is not None:
        # automatically generate the complementary sequence
        seq = []
        seq.append(args.seq_ds)

        comp = ''
        for i in reversed(range(len(seq[0]))):
            nuc = seq[0][i]
            if nuc == 'A':
                comp += 'T'
            elif nuc == 'T':
                comp += 'A'
            elif nuc == 'G':
                comp += 'C'
            elif nuc == 'C':
                comp += 'G'

        seq.append(comp)
        n_nt_all = 2 * len(args.seq_ds)
        n_seq = 2

    elif args.seqfile is not None:
        seq = ''
        for l in args.seqfile:
            if l[0] == '>' or l[0] == '#' or l[0] == ';':
                continue
            seq += l.strip()
        n_nt_all = len(seq)

    else:
        print ('Error: this error should be detected by mutually_exclusive_group of the parser')
        sys.exit(2)

    if args.circ:
        if args.end5 == 'S':
            print ('end5 option is ignored for circRNA')
            args.end5 = 'P'

    if args.end5 == 'P':
        offset_end5 = 1
    elif args.end5 == 'S':
        offset_end5 = 0
    else:
        print ('end5 option has to be Either S or P')
        sys.exit(2)

    ###########################
    ## Print the setting     ##
    ###########################
    print('#seq: {:d}'.format(n_seq))
    print('#nt: {:d}'.format(n_nt_all))
    print('Sequence:')
    for i in range(n_seq):
        print('    {:d} ({:d} nt.): {:s}'.format(i+1, len(seq[i]), seq[i]))

    if args.circ:
        print ('This is a circRNA. (--end5 option will be ignored)')
    else:
        print("This is a linear DNA with the 5'-end starting with %s" % (args.end5,))

    if args.exvfile:
        print('A file for exv pair list will be generated.')
    

    ################################################
    ## This is the main object to be constructed  ##
    ################################################
    ns = NinfoSet()
    
    ##############
    ## Bond     ##
    ##############

    nmp_pre = 0
    for iseq in range(n_seq):

        iunit = iseq + 1

        # for the first nt
        if args.end5 == 'P':
            # P - S
            bl = BondLength(iunit1=iunit, iunit2=iunit,
                            imp1=nmp_pre+0+offset_end5, imp2=nmp_pre+1+offset_end5, 
                            imp1un=0+offset_end5, imp2un=1+offset_end5,
                            native=BDNA.BL_PS, factor=1.0, correct_mgo=1.0, coef=CHT18.BL_PS, type_str='PS')
            ns.bondlengths.append(bl)
    
        # S - B
        native, coef, type_str = bl_native(iseq, 1)
        bl = BondLength(iunit1=iunit, iunit2=iunit,
                        imp1=nmp_pre+1+offset_end5, imp2=nmp_pre+2+offset_end5, 
                        imp1un=1+offset_end5, imp2un=2+offset_end5,
                        native=native, factor=1.0, correct_mgo=1.0, coef=coef, type_str=type_str)
        ns.bondlengths.append(bl)
        
        # for the second through last nt
        for i in range(2, len(seq[iseq])+1):
            imp_P = 3 * (i-1) + offset_end5
            imp_S0 = imp_P - 2
            imp_S = imp_P + 1
            imp_B = imp_S + 1
            # S0 - P (from previous nt)
            bl = BondLength(iunit1=iunit, iunit2=iunit,
                            imp1=nmp_pre+imp_S0, imp2=nmp_pre+imp_P, imp1un=imp_S0, imp2un=imp_P,
                            native=BDNA.BL_SP,factor=1.0,correct_mgo=1.0,coef=CHT18.BL_SP,type_str='SP')
            ns.bondlengths.append(bl)
        
            # P - S
            bl = BondLength(iunit1=iunit, iunit2=iunit,
                            imp1=nmp_pre+imp_P, imp2=nmp_pre+imp_S, imp1un=imp_P, imp2un=imp_S,
                            native=BDNA.BL_PS,factor=1.0,correct_mgo=1.0,coef=CHT18.BL_PS,type_str='PS')
            ns.bondlengths.append(bl)
        
            # S - B
            native, coef, type_str = bl_native(iseq, i)
            bl = BondLength(iunit1=iunit, iunit2=iunit,
                            imp1=nmp_pre+imp_S, imp2=nmp_pre+imp_B, imp1un=imp_S, imp2un=imp_B,
                            native=native, factor=1.0, correct_mgo=1.0, coef=coef, type_str=type_str)
            ns.bondlengths.append(bl)
    
        ## circRNA: connect the 5' and 3' ends
        if args.circ:
            imp_S0 = 3*len(seq[iseq]) - 1 # S of 3'-end
            imp_P  = 1          # New phosphate between 3'-end and 5'-end
    
            # S0 - P  (S of 3'-end to P of 5'-end)
            bl = BondLength(iunit1=1,iunit2=1,imp1=imp_S0,imp2=imp_P,imp1un=imp_S0,imp2un=imp_P,
                            native=BDNA.BL_SP,factor=1.0,correct_mgo=1.0,coef=CHT18.BL_SP,type_str='SP')
            ns.bondlengths.append(bl)

            nmp_pre += 3 * len(seq[iseq])
        
        else:
            if args.end5 == 'P':
                nmp_pre += 3 * len(seq[iseq])

            else:
                nmp_pre += 3 * len(seq[iseq]) - 1
    
    
    ##############
    ## Angle    ##
    ##############
    nmp_pre = 0
    for iseq in range(n_seq):

        iunit = iseq + 1
        n_nt = len(seq[iseq])

        for i in range(1, n_nt+1):
            imp_P  = 3*(i-1) + offset_end5   # This can be 0 if i == 1 and end5 == 'S'
            imp_S  = imp_P + 1
            imp_B  = imp_P + 2
            imp_P1 = imp_P + 3
            imp_S1 = imp_P + 4
    
            if i == n_nt and args.circ:
                imp_P1 = 1
                imp_S1 = 2
    
            """ P - S - B
            1st nucleotide: Only when there is P
            2nd through last: Always 
            """
            if (i > 1) or (i == 1 and args.end5 == 'P'):
                native, coef, type_str = ba_PSB_native(iseq, i)
                ba = BondAngle(iunit1=iunit, iunit2=iunit,
                               imp1=nmp_pre+imp_P, imp2=nmp_pre+imp_S, imp3=nmp_pre+imp_B, 
                               imp1un=imp_P, imp2un=imp_S, imp3un=imp_B, 
                               native=native, factor=1.0, correct_mgo=1.0, 
                               coef=coef, type_str=type_str)
                ns.bondangles.append(ba)
    
            """ P - S - P1
            1st nucleotide: Only when there is P
            2nd through 2nd last: Always 
            Last nucleotide: Only when circular
            """
            if (i > 1 and i < n_nt) or (i == 1 and args.end5 == 'P') or (i == n_nt and args.circ):
                ba = BondAngle(iunit1=iunit, iunit2=iunit,
                               imp1=nmp_pre+imp_P, imp2=nmp_pre+imp_S, imp3=nmp_pre+imp_P1, 
                               imp1un=imp_P, imp2un=imp_S, imp3un=imp_P1, 
                               native=BDNA.BA_PSP, factor=1.0, correct_mgo=1.0, 
                               coef=CHT18.BA_PSP, type_str='PSP')
                ns.bondangles.append(ba)
    
            """ B - S - P1, S - P1 - S1
            1st through 2nd last: Always 
            Last nucleotide: Only when circular
            """
            if (i < n_nt) or (i == n_nt and args.circ):
                native, coef, type_str = ba_BSP_native(iseq, i)
                ba = BondAngle(iunit1=iunit, iunit2=iunit, 
                               imp1=nmp_pre+imp_B, imp2=nmp_pre+imp_S, imp3=nmp_pre+imp_P1,
                               imp1un=imp_B, imp2un=imp_S, imp3un=imp_P1,
                               native=native,factor=1.0, correct_mgo=1.0,
                               coef=coef, type_str=type_str)
                ns.bondangles.append(ba)
    
                ba = BondAngle(iunit1=iunit, iunit2=iunit,
                               imp1=nmp_pre+imp_S, imp2=nmp_pre+imp_P1, imp3=nmp_pre+imp_S1,
                               imp1un=imp_S, imp2un=imp_P1, imp3un=imp_S1,
                               native=BDNA.BA_SPS, factor=1.0, correct_mgo=1.0,
                               coef=CHT18.BA_SPS, type_str='SPS')
                ns.bondangles.append(ba)

        if args.circ:
            nmp_pre += 3 * n_nt
        
        else:
            if args.end5 == 'P':
                nmp_pre += 3 * n_nt

            else:
                nmp_pre += 3 * n_nt - 1
    
    
    ##############
    ## Stack    ##
    ##############
    # for the second through one before last nt
    if args.end5 == 'P':
        nt_start = 2
    else:
        nt_start = 3

    nmp_pre = 0
    for iseq in range(n_seq):

        iunit = iseq + 1
        n_nt = len(seq[iseq])

        for i in range(nt_start, n_nt):
            imp_P2 = 3*(i-1) + offset_end5
            imp_P1 = imp_P2 - 3
            imp_S1 = imp_P2 - 2
            imp_B1 = imp_P2 - 1
            imp_S2 = imp_P2 + 1
            imp_B2 = imp_P2 + 2
            imp_P3 = imp_P2 + 3
        
            native, type_str = bs_native(iseq, i-1,i)
            h, s, Tm = CHT18.ST_U0[ type_str[0]+type_str[2] ]
            bs = BaseStackDT(iunit1=iunit, iunit2=iunit,
                               imp1=nmp_pre+imp_B1, imp2=nmp_pre+imp_B2,
                               imp1un=imp_B1, imp2un=imp_B2,
                               native=native, factor=0.0,correct_mgo=1.0,coef=CHT18.ST_DIST,type_str=type_str,
                               h=h, s=s, Tm=Tm,
                               dih1_imp1=nmp_pre+imp_P1, dih1_imp2=nmp_pre+imp_S1, dih1_imp3=nmp_pre+imp_P2, dih1_imp4=nmp_pre+imp_S2,
                               dih1_iunit1=iunit,dih1_iunit2=iunit,
                               dih1_imp1un=imp_P1, dih1_imp2un=imp_S1, dih1_imp3un=imp_P2, dih1_imp4un=imp_S2,
                               dih1_native=BDNA.DIH_PSPS,dih1_coef=CHT18.ST_DIH,dih1_type_str='PSPS',
                               dih2_imp1=nmp_pre+imp_S1, dih2_imp2=nmp_pre+imp_P2, dih2_imp3=nmp_pre+imp_S2, dih2_imp4=nmp_pre+imp_P3,
                               dih2_iunit1=iunit,dih2_iunit2=iunit,
                               dih2_imp1un=imp_S1, dih2_imp2un=imp_P2, dih2_imp3un=imp_S2, dih2_imp4un=imp_P3,
                               dih2_native=BDNA.DIH_SPSP,dih2_coef=CHT18.ST_DIH,dih2_type_str='SPSP')
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
            native, type_str = bs_native(iseq, n_nt-1,n_nt)
            bs = BaseStackDT(iunit1=iunit, iunit2=iunit,
                             imp1=nmp_pre+imp_B1, imp2=nmp_pre+imp_B2, imp1un=nmp_pre+imp_B1, imp2un=nmp_pre+imp_B2,
                             native=native, factor=0.0,correct_mgo=1.0,coef=CHT18.ST_DIST,type_str=type_str,
                             dih1_imp1=nmp_pre+imp_P1, dih1_imp2=nmp_pre+imp_S1, dih1_imp3=nmp_pre+imp_P2, dih1_imp4=nmp_pre+imp_S2,
                             dih1_iunit1=iunit, dih1_iunit2=iunit,
                             dih1_imp1un=imp_P1, dih1_imp2un=imp_S1, dih1_imp3un=imp_P2, dih1_imp4un=imp_S2,
                             dih1_native=BDNA.DIH_PSPS,dih1_coef=CHT18.ST_DIH,dih1_type_str='PSPS',
                             dih2_imp1=nmp_pre+imp_S1, dih2_imp2=nmp_pre+imp_P2, dih2_imp3=nmp_pre+imp_S2, dih2_imp4=nmp_pre+imp_P3,
                             dih2_iunit1=iunit, dih2_iunit2=iunit,
                             dih2_imp1un=imp_S1, dih2_imp2un=imp_P2, dih2_imp3un=imp_S2, dih2_imp4un=imp_P3,
                             dih2_native=BDNA.DIH_SPSP,dih2_coef=CHT18.ST_DIH,dih2_type_str='SPSP')
            ns.basestackDTs.append(bs)
    
            # (n_nt) : (1)
            imp_P1 = 3*n_nt - 2
            imp_S1 = 3*n_nt - 1
            imp_B1 = 3*n_nt
            imp_P2 = 1
            imp_S2 = 2
            imp_B2 = 3
            imp_P3 = 4
            native, type_str = bs_native(iseq, n_nt,1)
            bs = BaseStackDT(iunit1=iunit, iunit2=iunit,
                             imp1=nmp_pre+imp_B1, imp2=nmp_pre+imp_B2, imp1un=nmp_pre+imp_B1, imp2un=nmp_pre+imp_B2,
                             native=native, factor=0.0, correct_mgo=1.0, coef=CHT18.ST_DIST, type_str=type_str,
                             dih1_imp1=nmp_pre+imp_P1, dih1_imp2=nmp_pre+imp_S1, dih1_imp3=nmp_pre+imp_P2, dih1_imp4=nmp_pre+imp_S2,
                             dih1_iunit1=iunit,dih1_iunit2=iunit,
                             dih1_imp1un=imp_P1, dih1_imp2un=imp_S1, dih1_imp3un=imp_P2, dih1_imp4un=imp_S2,
                             dih1_native=BDNA.DIH_PSPS,dih1_coef=CHT18.ST_DIH,dih1_type_str='PSPS',
                             dih2_imp1=nmp_pre+imp_S1, dih2_imp2=nmp_pre+imp_P2, dih2_imp3=nmp_pre+imp_S2, dih2_imp4=nmp_pre+imp_P3,
                             dih2_iunit1=iunit,dih2_iunit2=iunit,
                             dih2_imp1un=imp_S1, dih2_imp2un=imp_P2, dih2_imp3un=imp_S2, dih2_imp4un=imp_P3,
                             dih2_native=BDNA.DIH_SPSP,dih2_coef=CHT18.ST_DIH,dih2_type_str='SPSP')
            ns.basestackDTs.append(bs)

        if args.circ:
            nmp_pre += 3 * n_nt
        
        else:
            if args.end5 == 'P':
                nmp_pre += 3 * n_nt

            else:
                nmp_pre += 3 * n_nt - 1

    
    ##############
    ## H-bond   ##
    ##############
    hblist = []
    if args.hbfile is not None:
        for l in args.hbfile:
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
    
    elif args.hb_canonical:

        for iseq in range(n_seq):
            for jseq in range(n_seq):
                for inuc_id, inuc in enumerate(seq[iseq]):
                    if inuc_id == 0 or inuc_id == len(seq[iseq])-1:
                        continue 

                    for jnuc_id, jnuc in enumerate(seq[jseq]):
                        if jnuc_id == 0 or jnuc_id == len(seq[jseq])-1:
                            continue 
                        
                        if iseq == jseq:
                            if inuc_id < jnuc_id + NNHB_NT_SEP:
                                continue
                        
                        if (inuc,jnuc) in (('A','T'), ('T','A')):
                            hblist.append(('CAN', inuc_id+1, 'B', jnuc_id+1, 'B', 2))
                            
                        elif (inuc,jnuc) in (('G','C'), ('C','G')):
                            hblist.append(('CAN', inuc_id+1, 'B', jnuc_id+1, 'B', 3))

        
    for c in hblist:
        nt_1 = c[1]
        mp_1 = c[2]
        nt_2 = c[3]
        mp_2 = c[4]
        nHB  = c[5]
    
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
             dih0_native, dih1_native, dih2_native, nHB) = hb_BDNA_native(0, nt_1, 0, nt_2)
    
#        elif c[0] == 'NON':  ## Other hydrogen bonds (including non-canonical basepairs)
#            xyz1 = chain.get_atom( imp_1 - 1 ).xyz
#            xyz2 = chain.get_atom( imp_2 - 1 ).xyz
#            xyz3 = chain.get_atom( imp_3 - 1 ).xyz
#            xyz4 = chain.get_atom( imp_4 - 1 ).xyz
#            xyz5 = chain.get_atom( imp_5 - 1 ).xyz
#            xyz6 = chain.get_atom( imp_6 - 1 ).xyz
#    
#            v12 = xyz1 - xyz2
#            v13 = xyz1 - xyz3
#            v53 = xyz5 - xyz3
#            v42 = xyz4 - xyz2
#            v46 = xyz4 - xyz6
#        
#            d1212 = v12.dot(v12)
#            dist_native = math.sqrt( d1212 )
#        
#            cos_theta = v13.dot(v12) / math.sqrt(v13.dot(v13) * d1212)
#            ang1_native = math.acos(cos_theta) / math.pi * 180.0
#        
#            cos_theta = v12.dot(v42) / math.sqrt(v42.dot(v42) * d1212)
#            ang2_native = math.acos(cos_theta) / math.pi * 180.0
#        
#            c4212 = v42.cross(v12)
#            c1213 = v12.cross(v13)
#            dih = math.atan2( v42.dot(c1213) * math.sqrt( v12.dot(v12)), c4212.dot(c1213))
#            dih0_native = dih / math.pi * 180.0
#    
#            m = v53.cross(v13)
#            n = c1213 * -1
#            dih = math.atan2( v53.dot(n) * math.sqrt( v13.dot(v13)), m.dot(n))
#            dih1_native = dih / math.pi * 180.0
#    
#            m = c4212 * -1
#            n = v42.cross(v46)
#            dih = math.atan2( v12.dot(n) * math.sqrt( v42.dot(v42)), m.dot(n))
#            dih2_native = dih / math.pi * 180.0
#    
        hb = HBondDT(iunit1=1,iunit2=1,
                     imp1=imp_1,imp2=imp_2,imp1un=imp_1,imp2un=imp_2,
                     native=dist_native, factor=CHT18.HB_U0*nHB, correct_mgo=1.0,coef=CHT18.HB_DIST,
    
                     ang1_imp1=imp_3, ang1_imp2=imp_1, ang1_imp3=imp_2, ang1_iunit1=1,ang1_iunit2=1,
                     ang1_imp1un=imp_3, ang1_imp2un=imp_1, ang1_imp3un=imp_2, 
                     ang1_native=ang1_native, ang1_coef=CHT18.HB_ANGL,
    
                     ang2_imp1=imp_4, ang2_imp2=imp_2, ang2_imp3=imp_1, ang2_iunit1=1,ang2_iunit2=1,
                     ang2_imp1un=imp_4, ang2_imp2un=imp_2, ang2_imp3un=imp_1, 
                     ang2_native=ang2_native, ang2_coef=CHT18.HB_ANGL,
    
                     dih0_imp1=imp_3, dih0_imp2=imp_1, dih0_imp3=imp_2, dih0_imp4=imp_4,dih0_iunit1=1,dih0_iunit2=1,
                     dih0_imp1un=imp_3, dih0_imp2un=imp_1, dih0_imp3un=imp_2, dih0_imp4un=imp_4,
                     dih0_native=dih0_native,dih0_coef=CHT18.HB_DIH_HBOND,
    
                     dih1_imp1=imp_2, dih1_imp2=imp_1, dih1_imp3=imp_3, dih1_imp4=imp_5,dih1_iunit1=1,dih1_iunit2=1,
                     dih1_imp1un=imp_2, dih1_imp2un=imp_1, dih1_imp3un=imp_3, dih1_imp4un=imp_5,
                     dih1_native=dih1_native,dih1_coef=CHT18.HB_DIH_CHAIN,
    
                     dih2_imp1=imp_1, dih2_imp2=imp_2, dih2_imp3=imp_4, dih2_imp4=imp_6,dih2_iunit1=1,dih2_iunit2=1,
                     dih2_imp1un=imp_1, dih2_imp2un=imp_2, dih2_imp3un=imp_4, dih2_imp4un=imp_6,
                     dih2_native=dih2_native,dih2_coef=CHT18.HB_DIH_CHAIN,
                     )
    
        ns.hbondDTs.append(hb)
    

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
    if args.exvfile:
        n_sep_nlocal_P = 3
        n_sep_nlocal_S = 3
        n_sep_nlocal_B = 2

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
                    imp_plus_sep[imp] = imp + n_sep_nlocal_B
                elif itype == 1: # P
                    imp_plus_sep[imp] = imp + n_sep_nlocal_P
                else: # S
                    imp_plus_sep[imp] = imp + n_sep_nlocal_S
            else:
                if itype == 0: # P
                    imp_plus_sep[imp] = imp + n_sep_nlocal_P
                elif itype == 1: # S
                    imp_plus_sep[imp] = imp + n_sep_nlocal_S
                else: # B
                    imp_plus_sep[imp] = imp + n_sep_nlocal_B

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

