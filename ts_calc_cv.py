#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2014/06/17
@author: Naoto Hori
'''

import sys
from cafysis.file_io.ts import TsFile
import glob

BOLTZC = 0.0019872041
#E_VELO_CUTOFF = 160.0
#E_CUTOFF = 110.0

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print 'Usage: %SCRIPT [search query (use "")] [skip] [end] [Energy cutoff]'
        sys.exit(2)

    nframe_skip = int(sys.argv[2])
    nframe_end = int(sys.argv[3])
    E_CUTOFF = float(sys.argv[4])

    for tsfile in glob.glob(sys.argv[1]):
        ts = TsFile(tsfile)
        ts.open_to_read()
        ts.read_header()
        
        e_sq = 0.0
        e = 0.0
        n = 0
        first = True
        temp = 0.0
        iframe = 0
        while ts.has_more_data():

            tsdata, tslines = ts.read_onestep()
            iframe += 1
            if iframe <= nframe_skip:
                continue
            if iframe > nframe_end:
                break

            if first:
                temp = float(tsdata[0][ts.head_col.temp])
                first = False
            else:
                if temp != float(tsdata[0][ts.head_col.temp]):
                    print ('temperature is not identical')
                    sys.exit(2)

            etot = float(tsdata[0][ts.head_col.e_tot])
            #evelo = float(tsdata[0][ts.head_col.e_velo])
            if etot > E_CUTOFF:
                continue

            e_sq += etot * etot
            e += etot
            n += 1

        fn = float(n)
        print temp, ((e_sq/fn - (e/fn)**2) / (BOLTZC*temp*temp))

