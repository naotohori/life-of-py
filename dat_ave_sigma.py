#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2014/02/17
@author: Naoto Hori
'''

import sys
import glob
import math

COLUM_DIST = 4 - 1

if __name__ == '__main__':
    if not len(sys.argv) in (4,5):
        print('Usage: % SCRIPT [GLOB query ("" needed)] [#step to be ignored] [output]')
        print('  or : % SCRIPT [GLOB query ("" needed)] [#step to be ignored] [final step] [output]')
        sys.exit(2)

    files = glob.glob(sys.argv[1])
    step_ignore = int(sys.argv[2])
    flg_final = False
    if len(sys.argv) == 5:
        flg_final = True
        step_final  = int(sys.argv[3])
    filepath_out = sys.argv[-1]

    f_out = open(filepath_out,'w')
    for filepath in files:
        f_in = open(filepath,'r')

        step = 0
        ave = 0.0
        ave2 = 0.0
        num = 0
        for line in f_in:
            step = step + 1
            if step < step_ignore:
                continue 

            linesp = line.split()
            dist = float(linesp[COLUM_DIST])
            ave = ave + dist
            ave2 = ave2 + dist * dist
            num = num + 1

            if flg_final and step == step_final:
                break

        f_in.close()

        if num != (step - step_ignore +1):
            print(('Warning: number step is %i in %s' % (num, filepath)))
        if num == 0:
            continue
            
        ave = ave / num
        ave2 = ave2 / num
        sigma = math.sqrt(ave2 - ave * ave)

        ## filepath is like "/Users/hori/bwyv/cM500/500_00100_002/Energies.dat"
        ## Here "500", "00100", and "002" should be extracteed.
        cM    = filepath.replace('/','_').split('_')[-4]
        force = filepath.replace('/','_').split('_')[-3]
        runid = filepath.replace('/','_').split('_')[-2]

        f_out.write('%s %s %s %6.2f %6.2f\n' % (cM, force, runid, ave, sigma))

    f_out.close()


