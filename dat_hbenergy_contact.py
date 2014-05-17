#!/usr/bin/env python

# IDX TYP  CG1  NT1   CG2  NT2   DST    REF
#   1 CAN    11 C03B    56 G18B   5.746   5.655
#CLM_HB_IDX = 1 - 1
#CLM_HB_TYP = 2 - 1
#CLM_HB_CG1 = 3 - 1
#CLM_HB_NT1 = 4 - 1
#CLM_HB_CG2 = 5 - 1
#CLM_HB_NT2 = 6 - 1
#CLM_HB_DST = 7 - 1
#CLM_HB_REF = 8 - 1
#CLM_HB_2ND = 9 - 1
CLM_HB_NBOND = 10 - 1

#CUTOFF_CONTACT_LOW = 0.8
#CUTOFF_CONTACT_HIG = 1.2

import sys
import os
import glob

if len(sys.argv) != 4:
    print ('Usage: SCRIPT [HB file (bwyv_forEne2.hb)] [dir_search] [CUTOFF_CONTACT_ENE]')
    sys.exit(2)

filepath_hb = sys.argv[1]
dir_search = sys.argv[2]
CUTOFF_CONTACT_ENE = float(sys.argv[3])

filename_out = 'hbcon2_%.2f.out' % (CUTOFF_CONTACT_ENE,)

f_hb = open(filepath_hb, 'r')
hbs = []
for l in f_hb:
    hbs.append(l.split())
        
''' Collect dirnames and detect Ion conc. and Forces. '''
#dirs = glob.glob('*/*_*_*/')
if dir_search[-1] != '/':
    dir_search = dir_search + '/'
dirs = glob.glob(dir_search)
simulations = []
for d in dirs:
    d_sp = d.split('/')[-2].split('_')
    cM  = d_sp[-3]
    frc = d_sp[-2]
    rnd = d_sp[-1]
    type_pNcM = d.split('/')[-3][0:2]
    if type_pNcM not in ('pN','cM'):
        print 'Error: not pN nor cM'
        sys.exit(2)
    simulations.append((type_pNcM,cM,frc,rnd))


''' Loop for each simulation set. '''
orig_dir = os.getcwd()
for sim in simulations:
    type_pNcM, cM, frc, rnd = sim
    os.chdir(orig_dir)
    os.chdir('%s/%s%s/%s_%s_%s' % (orig_dir, type_pNcM, cM, cM, frc, rnd))

    if not os.path.exists('Energies2_HB.dat'):
        print ('Warning: No file Energies2_HB.dat in (type_pNcM, cM, frc, rnd)=(%s,%s,%s,%s)' 
              % (type_pNcM, cM, frc, rnd) )
        continue

    if os.path.exists(filename_out):
        if os.stat(filename_out).st_mtime >= os.stat('Energies2_HB.dat').st_mtime:
            continue

    ''' Open output file'''
    f_out = open(filename_out, 'w')

    for l in open('Energies2_HB.dat'):

        l = l.split()
        if len(l) != len(hbs):
            print ('Error: len(l) != len(hbs)')
            sys.exit(2)

        ''' Prepare a dictionary for counting'''
        contact_seq = []
    
        ''' Loop for HB pairs.'''
        for ihb, hb in enumerate(hbs):
            #idx = int(hb[CLM_HB_IDX])
            #cg1 = int(hb[CLM_HB_CG1])
            #cg2 = int(hb[CLM_HB_CG2])
            #nt1 = hb[CLM_HB_NT1]
            #nt2 = hb[CLM_HB_NT2]
            #ref = float(hb[CLM_HB_REF])
            #name2nd = hb[CLM_HB_2ND]  # name of the 2ndary structure
            nbond = int(hb[CLM_HB_NBOND])
    
            ene = float(l[ihb])
            if ene <= CUTOFF_CONTACT_ENE * nbond:
                contact_seq.append('1')
            else:
                contact_seq.append('0')
    
        for con in contact_seq:
            f_out.write('%1s' % (con,))
        f_out.write('\n')
