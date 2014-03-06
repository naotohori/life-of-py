#!/usr/bin/env python

# IDX TYP  CG1  NT1   CG2  NT2   DST    REF
#   1 CAN    11 C03B    56 G18B   5.746   5.655
#CLM_HB_IDX = 1 - 1
#CLM_HB_TYP = 2 - 1
#CLM_HB_CG1 = 3 - 1
CLM_HB_NT1 = 4 - 1
#CLM_HB_CG2 = 5 - 1
CLM_HB_NT2 = 6 - 1
#CLM_HB_DST = 7 - 1
CLM_HB_REF = 8 - 1
CLM_HB_2ND = 9 - 1

CUTOFF_CONTACT = 1.2

import sys
import os
import glob

if len(sys.argv) != 2:
    print 'Usage: SCRIPT [HB file (bwyv.hb)]' 
    sys.exit(2)

filepath_hb = sys.argv[1]

f_hb = open(filepath_hb, 'r')
hbs = []
for l in f_hb:
    hbs.append(l.split())
        
''' Collect dirnames and detect Ion conc. and Forces. '''
dirs = glob.glob('*/*_*_*/')
simulations = []
for d in dirs:
    d_sp = d.split('/')[-2].split('_')
    cM  = d_sp[-3]
    frc = d_sp[-2]
    rnd = d_sp[-1]
    simulations.append((cM,frc,rnd))


''' Loop for each simulation set. '''
orig_dir = os.getcwd()
for sim in simulations:
    cM, frc, rnd = sim
    os.chdir('%s/cM%s/%s_%s_%s' % (orig_dir, cM, cM, frc, rnd))

    ''' Open distance files'''
    files_data = []
    exist_all_files = True
    for hb in hbs:
        nt1 = hb[CLM_HB_NT1]
        nt2 = hb[CLM_HB_NT2]
        datafile = 'dist_%s_%s.out' % (nt1,nt2)
        try:
            files_data.append(open(datafile,'r'))
        except:
            exist_all_files = False
            continue

    if not exist_all_files:
        continue

    ''' Count number of lines'''
    num_lines = 0
    for line in files_data[0]:
        if line.find('#') != -1:
            continue
        num_lines = num_lines + 1
    files_data[0].seek(0)

    ''' Open output file'''
    f_out = open('hb2nd.out', 'w')

    iline = 0
    while (iline < num_lines):
        iline = iline + 1

        ''' Prepare a dictionary for counting'''
        contact2nd = {'S1':0, 'S2':0, 'L1':0, 'L2':0}
    
        ''' Loop for HB pairs.'''
        for ihb, hb in enumerate(hbs):
            #idx = int(hb[CLM_HB_IDX])
            #cg1 = int(hb[CLM_HB_CG1])
            #cg2 = int(hb[CLM_HB_CG2])
            #nt1 = hb[CLM_HB_NT1]
            #nt2 = hb[CLM_HB_NT2]
            ref = float(hb[CLM_HB_REF])
            name2nd = hb[CLM_HB_2ND]  # name of the 2ndary structure
            contact = CUTOFF_CONTACT * ref
    
            ''' Calc statistical data (percentage).'''
            line = files_data[ihb].readline()
            while line.find('#') != -1:
                line = files_data[ihb].readline()
    
            d = float(line.strip())
            if d < contact:
                contact2nd[ name2nd ] = contact2nd[ name2nd ] + 1         
    
        f_out.write('%3i %3i %3i %3i\n' % (contact2nd['S1'],
                                           contact2nd['S2'],
                                           contact2nd['L1'],
                                           contact2nd['L2']) )
    
    os.chdir(orig_dir)

