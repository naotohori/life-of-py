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
CLM_HB_2ND = 9 - 1

import sys
import os
import glob

if len(sys.argv) != 5:
    print(('Usage: SCRIPT [HB file (bwyv_forEne2.hb)] [dir_search]'
                        +' [filename con (hbcon2_-1.47.out)]'
                        +' [filename hb2nd (hb2nd_-1.47.out)]'))
    sys.exit(2)

filepath_hb = sys.argv[1]
dir_search = sys.argv[2]
filename_con = sys.argv[3]
filename_hb2nd = sys.argv[4]

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
    type = d.split('/')[-3][0:2]
    if type not in ('pN','cM'):
        print('Error: not pN nor cM')
        sys.exit(2)
    simulations.append((type,cM,frc,rnd))


''' Loop for each simulation set. '''
orig_dir = os.getcwd()
for sim in simulations:
    type, cM, frc, rnd = sim
    os.chdir(orig_dir)
    os.chdir('%s/%s%s/%s_%s_%s' % (orig_dir, type, cM, cM, frc, rnd))

    if not os.path.exists(filename_con):
        continue

    if os.path.exists(filename_hb2nd):
        if os.stat(filename_hb2nd).st_mtime >= os.stat(filename_con).st_mtime:
            continue

    ''' Open output file'''
    f_out = open(filename_hb2nd, 'w')

    for l in open(filename_con,'r'):

        ''' Prepare a dictionary for counting'''
        contact2nd = {'S1':0, 'S2':0, 'L1':0, 'L2':0}

        l = l.strip()
        if len(l) != len(hbs):
            print ('Error: len(l) != len(hbs)')
            sys.exit(2)
    
        ''' Loop for HB pairs.'''
        for ihb, hb in enumerate(hbs):
    
            if l[ihb] == '1':
                name2nd = hb[CLM_HB_2ND]  # name of the 2ndary structure
                contact2nd[name2nd] += 1
    
        f_out.write('%3i %3i %3i %3i\n' % (contact2nd['S1'],
                                           contact2nd['S2'],
                                           contact2nd['L1'],
                                           contact2nd['L2']) )
    
os.chdir(orig_dir)
