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

CUTOFF_CONTACT_LOW = 0.8
CUTOFF_CONTACT_HIG = 1.2

import sys
import os
import glob

if len(sys.argv) == 5:
    step_final  = int(sys.argv[-1])
    flg_final = True
elif len(sys.argv) == 4:
    flg_final = False
else:
    print('Usage: SCRIPT [HB file (bwyv.hb)] [dir_search] [step_ignore]')
    print(' or  : SCRIPT [HB file (bwyv.hb)] [dir_search] [step_ignore] [step_final]')
    sys.exit(2)

filepath_hb = sys.argv[1]
dir_search = sys.argv[2]
step_ignore = int(sys.argv[3])


f_hb = open(filepath_hb, 'r')
hbs = []
idx = []
cg1 = []
cg2 = []
nt1 = []
nt2 = []
ref = []
num_con = 0
for l in f_hb:
    l = l.split()
    idx.append(int(l[CLM_HB_IDX]))
    cg1.append(int(l[CLM_HB_CG1]))
    cg2.append(int(l[CLM_HB_CG2]))
    nt1.append(l[CLM_HB_NT1])
    nt2.append(l[CLM_HB_NT2])
    ref.append(float(l[CLM_HB_REF]))
    num_con = num_con + 1
f_hb.close()
        
''' Collect dirnames and detect Ion conc. and Forces. '''
dirs = glob.glob(dir_search)
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
    os.chdir('cM%s/%s_%s_%s' % (cM, cM, frc, rnd))

    f_data = []
    flg_skip = False
    for i in range(num_con):
        filename = 'dist_%s_%s.out' % (nt1[i], nt2[i])
        try:
            f_data.append(open(filename, 'r'))
        except:
            flg_skip = True
            continue

    if flg_skip:
        print(('Skip %s' % (os.getcwd(),) ))
        os.chdir(orig_dir)
        continue

    con = []
    for _ in range(num_con):
        con.append([])

    step = 0
    for l in f_data[0]:
        lines = [l]
        for i in range(1,num_con):
            lines.append(next(f_data[i]))

        if (lines[0].find('#') != -1):
            continue

        step = step + 1
        if step < step_ignore:
            continue

        for i in range(num_con):
            if CUTOFF_CONTACT_LOW * ref[i] <= float(lines[i].strip()) <= CUTOFF_CONTACT_HIG * ref[i]:
                con[i].append(1.0)
            else:
                con[i].append(0.0)

        if flg_final and step == step_final:
            break

    n = len(con[0])
    avg = []
    for c in con:
        avg.append(sum(c)/float(n))

    var = []
    for i,c in enumerate(con):
        var.append( sum( [pow(x - avg[i], 2) for x in c] ) / n )
    
    std = []
    for v in var:
        std.append( pow(v, 0.5) )

    cov = []
    for _ in range(num_con):
        cov.append([0.0] * num_con)

    for i in range(num_con):
        for j in range(i+1, num_con):
            covsum = 0.0
            for ii in range(n):
                covsum += (con[i][ii] - avg[i]) * (con[j][ii] - avg[j])
            cov[i][j] = covsum / float(n)

    cor = []
    for _ in range(num_con):
        cor.append([0.0] * num_con)

    for i in range(num_con):
        if std[i] == 0.0:
            continue
        for j in range(i+1, num_con):
            if std[j] == 0.0:
                continue
            cor[i][j] =  cov[i][j] / (std[i] * std[j])

    f_out = open('hb_cor.out', 'w')
    for i in range(num_con):
        for j in range(num_con):
            f_out.write('%i %i %f\n' % (i+1,j+1,cor[i][j]))
        f_out.write('\n\n')

    f_out.close()
    os.chdir(orig_dir)

