#!/usr/bin/env python

import math

cM_values = []
frc_values = []
#ihb_values = range(1,24)

#count = {'F':0, 'S1S2':0, 'S1S2L2':0, 'S1S2':0, 'S1':0, 'S2':0, 'U':0, 'other':0} 
states = {0:'F', 1:'I1', 2:'I2', 3:'U', 4:'other'}
sum_ratio = {}
sum2_ratio = {}
n_ratio = {}
sim_done = set()

import sys
if len(sys.argv) != 3:
    print('Usage: SCRIPT [input file (hb2nd.out)] [output dir (hb2nd/)]')
    sys.exit(2)

filepath_inp = sys.argv[1]
dirpath_out = sys.argv[2]

for l in open(filepath_inp,'r'):
    if l.find('#') != -1:
        continue

    l = l.split()

    cM = float(l[0])
    frc = float(l[1])
    rnd = l[2]

    sim = (cM, frc)
    sim_done.add(sim)

    if sim in n_ratio:
        n_ratio[sim] = n_ratio[sim] + 1
        for i in range(5):  
            x = float(l[i+3])
            sum_ratio[sim][i] += x
            sum2_ratio[sim][i] += x*x
    else:
        n_ratio[sim] = 1
        sum_ratio[sim] = []
        sum2_ratio[sim] = []
        for i in range(5):  
            x = float(l[i+3])
            sum_ratio[sim].append(x)
            sum2_ratio[sim].append(x*x)

    if cM not in cM_values:
        cM_values.append(cM)
    if frc not in frc_values:
        frc_values.append(frc)

cM_values.sort()
frc_values.sort()

for cM in cM_values:
    for frc in frc_values:
        sim = (cM,frc)
        if sim not in sim_done:
            continue
            
        f_out = open('%s/%04i_%02i.hb2nd' % (dirpath_out, cM, frc), 'w')
        f_out.write('#%i\n' % (n_ratio[sim],))
        for i in range(5):
            ave = sum_ratio[sim][i] / n_ratio[sim]
            #sd = math.sqrt( sum2_ratio[sim][i] / n_ratio[sim]  - ave * ave )
            se = math.sqrt( (sum2_ratio[sim][i] / float(n_ratio[sim]) - ave * ave )
                           / float(n_ratio[sim]) )
            f_out.write('%6s %7.3f %7.3f\n' % (states[i], ave, se ))

        f_out.close()
