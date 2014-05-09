#!/usr/bin/env python

# IDX TYP  CG1  NT1   CG2  NT2   DST    REF
#   1 CAN    11 C03B    56 G18B   5.746   5.655
CLM_HB_IDX = 1 - 1
CLM_HB_TYP = 2 - 1
CLM_HB_CG1 = 3 - 1
CLM_HB_NT1 = 4 - 1
CLM_HB_CG2 = 5 - 1
CLM_HB_NT2 = 6 - 1
CLM_HB_DST = 7 - 1
CLM_HB_REF = 8 - 1
CLM_HB_2ND = 9 - 1

CUTOFF_CONTACT = 1.2

import sys
import os
import glob

if len(sys.argv) == 5:
    step_final  = int(sys.argv[-1])
    flg_final = True
elif len(sys.argv) == 4:
    flg_final = False
else:
    print 'Usage: SCRIPT [HB file (bwyv.hb)] [output file] [step_ignore]'
    print ' or  : SCRIPT [HB file (bwyv.hb)] [output file] [step_ignore] [step_final]'
    sys.exit(2)

filepath_hb = sys.argv[1]
filepath_out = sys.argv[2]
step_ignore = int(sys.argv[3])


f_hb = open(filepath_hb, 'r')
hbs = []
for l in f_hb:
    hbs.append(l.split())
        
''' Collect dirnames and detect Ion conc. and Forces. '''
dirs = glob.glob('cM*/*_*_*/')
simulations = []
for d in dirs:
    d_sp = d.split('/')[-2].split('_')
    cM  = d_sp[-3]
    frc = d_sp[-2]
    rnd = d_sp[-1]
    simulations.append((cM,frc,rnd))

f_out = open(filepath_out, 'w')

''' Loop for each simulation set. '''
orig_dir = os.getcwd()
for sim in simulations:
    cM, frc, rnd = sim
    os.chdir('cM%s/%s_%s_%s' % (cM, cM, frc, rnd))

    ''' Loop for HB pairs.'''
    for hb in hbs:
        idx = int(hb[CLM_HB_IDX])
        cg1 = int(hb[CLM_HB_CG1])
        cg2 = int(hb[CLM_HB_CG2])
        nt1 = hb[CLM_HB_NT1]
        nt2 = hb[CLM_HB_NT2]
        ref = float(hb[CLM_HB_REF])
        #2nd = hb[CLM_HB_2ND]  # name of the 2ndary structure
        contact = CUTOFF_CONTACT * ref

        datafile = 'dist_%s_%s.out' % (nt1,nt2)
        try:
            f_data = open(datafile,'r')
        except:
            continue

        ''' Calc statistical data (percentage).'''
        step = 0
        num = 0 
        n_con = 0
        for line in f_data:
            if (line.find('#') != -1):
                continue

            step = step + 1
            if step < step_ignore:
                continue

            num = num + 1
            d = float(line.strip())
            if d < contact:
                n_con = n_con + 1

            if flg_final and step == step_final:
                break

        if num == 0:
            print ('Warning: %s has no data.' % (datafile,))
        else:
            ratio = 100.0 * float(n_con) / float(num)
            f_out.write('%s %s %s %5i %6.3f\n' % (cM,frc,rnd,idx,ratio))

    os.chdir(orig_dir)

