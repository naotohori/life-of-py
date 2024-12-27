#!/usr/bin/env python
'''
Created on 2024/10/06
@author: Huong
To get F-x paramenters if forgot to add twz into the input.toml file.
#  (1)step  (2)id (3)d[traps]   (4)d[beads]   (5)d[t1-b1]   (6)d[t2-b2]   (7)f[b1]      (8)f[b2]
        10000   1   15.305719       13.673979       3.0965533       1.8855414       64.541725       39.300501
        20000   1   15.320719       16.975788      0.67501561       2.7001408       14.069408       56.279267
        30000   1   15.335719       14.521236       2.9823167       2.0597383       62.160683       42.931302
        40000   1   15.350719       11.650796       2.5086845       3.1930296       52.288726       66.552589
        50000   1   15.365719       13.680551       2.4918765       1.5111290       51.938396       31.496591

'''

import sys
import os
import math
from lop.file_io.dcd import DcdFile

if len(sys.argv) != 4:
    print(' Usage: % SCRIPT [input DCD] [out log file] [output PDB] ')
    sys.exit(2)
dcd_file = sys.argv[1]
dcd = DcdFile(dcd_file)

out_file = sys.argv[2]
with open(out_file, 'r') as f:
    for l in f:
        if "## pair, imp1, imp2: " in l:
            print(l)
            pair = int(l.strip().split()[-3])
            id1 = int(l.strip().split()[-2]) - 1
            id2 = int(l.strip().split()[-1]) - 1
            print(pair)
            print(id1)
            print(id2)
        elif "## initial position 1: " in l:
            print(l)
            ini_t1=[float(x) for x in l.strip().split()[-3:]]
            print(ini_t1)
        elif "## initial position 2: " in l:
            print(l)
            ini_t2=[float(x) for x in l.strip().split()[-3:]]
            print(ini_t2)
        elif "## velocity 1: " in l:
            print(l)
            v_t1=[float(x) for x in l.strip().split()[-3:]]
            print(str(v_t1) + ' in A/step')
        elif "## velocity 2: " in l:
            print(l)
            v_t2=[float(x) for x in l.strip().split()[-3:]]
            print(str(v_t2) + ' in A/step')
        elif "Step" in l:
            print(l)
            int_step = int(l.strip().split()[-1])
            print(int_step)
        elif "# MD nstep_save: " in l:
            print(l)
            frame_freq = int(l.strip().split()[-1])
            print(frame_freq)
        # elif "Spring constant" current format did not write out the spring constant. So will have to get that from the name of the file now


print('spring constant')
k = float(dcd_file.split('_')[1])
print(str(k)+' kcal/mol/A^2 ')

k = k/0.0003*0.2/10
print(str(k)+' pN/A')

dcd.open_to_read()
dcd.read_header()

f_out = open(sys.argv[-1],'w')
f_out.write("#  (1)step  (2)id (3)d[traps]   (4)d[beads]   (5)d[t1-b1]   (6)d[t2-b2]   (7)f[b1]      (8)f[b2]\n")
nframe = 0
while dcd.has_more_data() :
    data = dcd.read_onestep()
    nframe += 1
    # (1)step
    s = int_step + nframe * frame_freq

    # (2)id is the same with pair
    # (3)d[traps]
    x_t1 = [ini + v*nframe*frame_freq for ini, v in zip(ini_t1, v_t1) ]
    x_t2 = [ini + v*nframe*frame_freq for ini, v in zip(ini_t2, v_t2) ]

    dt = math.sqrt( (x_t1[0] - x_t2[0]) ** 2
                  +(x_t1[1] - x_t2[1]) ** 2
                  +(x_t1[2] - x_t2[2]) ** 2 )

    # end-to-end distance which is (4)d[beads]
    d = math.sqrt( (data[id1][0] - data[id2][0]) ** 2
                  +(data[id1][1] - data[id2][1]) ** 2
                  +(data[id1][2] - data[id2][2]) ** 2 )

    # (5)d[t1-b1]
    d1 = math.sqrt( (x_t1[0] - data[id1][0]) **2
                   +(x_t1[1] - data[id1][1]) **2
                   +(x_t1[2] - data[id1][2]) **2 )

    # (6)d[t2-b2]
    d2 = math.sqrt( (x_t2[0] - data[id2][0]) **2
                   +(x_t2[1] - data[id2][1]) **2
                   +(x_t2[2] - data[id2][2]) **2 )      

    # (7)f[b1]
    f1 = k * d1 

    # (8)f[b2]
    f2 = k * d2   

    f_out.write('%15d %3d %20.6f  %20.6f %20.6f  %20.6f %20.6f  %20.6f \n' % (s, pair, dt, d, d1, d2, f1, f2))
dcd.close()

f_out.close()
