#!/usr/bin/env python

import math
import glob
import sys 
import os

S1_FORMED = 4
S2_FORMED = 4
L1L2_FORMED = 10

if len(sys.argv) == 4:
    step_ignore = int(sys.argv[2])
    step_final  = int(sys.argv[-1])
    flg_final = True
elif len(sys.argv) == 3:
    step_ignore = int(sys.argv[2])
    flg_final = False
else:
    print 'Usage: SCRIPT [output file] [step_ignore]'
    print ' or  : SCRIPT [output file] [step_ignore] [step_final]'
    sys.exit(2)

f_out = open(sys.argv[1], 'w')
f_out.write('#1  2  3     4       5       6       7      8\n')
f_out.write('#cM fr rn  Folded    S1S2    S1   Unfold   other\n')

''' Collect dirnames and detect Ion conc. and Forces. '''
dirs = glob.glob('cM*/*_*_*/')
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

    try:
        f = open('hb2nd.out','r')
    except:
        continue

    count = {'F':0, 'S1S2L1':0, 'S1S2L2':0, 'S1S2':0, 'S1':0, 'S2':0, 'U':0, 'other':0} 
    num = 0
    step = 0
    for l in f:
        if l.find('#') != -1:
            continue
        step = step + 1
        if step < step_ignore:
            continue

        l = l.split()
        S1 = False
        S2 = False
        L1L2 = False
        if int(l[0]) >= S1_FORMED:
            S1 = True
        if int(l[1]) >= S2_FORMED:
            S2 = True
        if (int(l[2]) + int(l[3])) >= L1L2_FORMED:
            L1L2 = True

        if S1 and S2 and L1L2:
            count['F'] = count['F'] + 1
        elif S1 and S2:
            count['S1S2'] = count['S1S2'] + 1
        elif S1:
            count['S1'] = count['S1'] + 1
        #elif S2:
        #    count['S2'] = count['S2'] + 1
        elif (not S1) and (not S2) and (not L1L2):
            count['U'] = count['U'] + 1
        else:
            count['other'] = count['other'] + 1
            #print 'Error hendesuyo'

        num = num + 1
        if flg_final and step == step_final:
            break

    f_out.write('%3s %2s %2s %7.3f %7.3f %7.3f %7.3f %7.3f\n' %
               (cM, frc, rnd, 100*count['F']      / float(num),
                              100*count['S1S2']   / float(num),
                              100*count['S1']     / float(num),
                              #100*count['S2']     / float(num),
                              100*count['U']      / float(num),
                              100*count['other']  / float(num))  )

f_out.close()

