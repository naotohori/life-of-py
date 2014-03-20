#!/usr/bin/env python

import math

cM_values = []
frc_values = []
#ihb_values = range(1,24)

n_ratio = {}
sum_ratio = {}
sum2_ratio = {}
sim_done = set()

import sys
if len(sys.argv) != 3:
    print ('Usage: SCRIPT [input file (hbprob.out)] [output dir (hbprob/)]')
    sys.exit(2)

filepath_inp = sys.argv[1]
dirpath_out = sys.argv[2]

for l in open(filepath_inp, 'r'):
    l = l.split()

    cM = float(l[0])
    frc = float(l[1])
    rnd = l[2]
    ihb = int(l[3])
    ratio = float(l[4])

    sim = (cM, frc)
    key = (sim, ihb)

    sim_done.add(sim)

    if n_ratio.has_key(key):
        n_ratio[key] = n_ratio[key] + 1
        sum_ratio[key] = sum_ratio[key] + ratio
        sum2_ratio[key] = sum2_ratio[key] + ratio * ratio
    else:
        n_ratio[key] = 1
        sum_ratio[key] = ratio
        sum2_ratio[key] = ratio * ratio

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
            
        f_out = open('%s/%04i_%02i.hb_prob' % (dirpath_out, cM, frc),'w')
        # Stem 1
        for ihb in range(1,6):
            key = (sim, ihb)
            ave = sum_ratio[key] / n_ratio[key]
            sd = math.sqrt( sum2_ratio[key] / n_ratio[key]  - ave * ave )
            f_out.write('%3i %7.3f %7.3f %3i\n' % (ihb, ave, sd, n_ratio[key]))
        f_out.write('\n\n')
        # Stem 2
        for ihb in range(6,9):
            key = (sim, ihb)
            ave = sum_ratio[key] / n_ratio[key]
            sd = math.sqrt( sum2_ratio[key] / n_ratio[key]  - ave * ave )
            f_out.write('%3i %7.3f %7.3f %3i\n' % (ihb, ave, sd, n_ratio[key]))
        f_out.write('\n\n')
        # C14-A25
        for ihb in range(9,11):
            key = (sim, ihb)
            ave = sum_ratio[key] / n_ratio[key]
            sd = math.sqrt( sum2_ratio[key] / n_ratio[key]  - ave * ave )
            f_out.write('%3i %7.3f %7.3f %3i\n' % (ihb, ave, sd, n_ratio[key]))
        f_out.write('\n\n')
        # L2
        for ihb in range(11,21):
            key = (sim, ihb)
            ave = sum_ratio[key] / n_ratio[key]
            sd = math.sqrt( sum2_ratio[key] / n_ratio[key]  - ave * ave )
            f_out.write('%3i %7.3f %7.3f %3i\n' % (ihb, ave, sd, n_ratio[key]))
        f_out.write('\n\n')
        # L1
        for ihb in range(21,24):
            key = (sim, ihb)
            ave = sum_ratio[key] / n_ratio[key]
            sd = math.sqrt( sum2_ratio[key] / n_ratio[key]  - ave * ave )
            f_out.write('%3i %7.3f %7.3f %3i\n' % (ihb, ave, sd, n_ratio[key]))

        f_out.close()
