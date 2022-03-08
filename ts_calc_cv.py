#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2014/06/17
@author: Naoto Hori
'''

import sys
from lop.file_io.ts import TsFile
import glob

BOLTZC = 0.0019872041
#E_VELO_CUTOFF = 160.0
#E_CUTOFF = 110.0

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print('Usage: %SCRIPT [search query (use "")] [skip] [end] [Energy cutoff]')
        sys.exit(2)

    nframe_skip = int(sys.argv[2])
    nframe_end = int(sys.argv[3])
    E_CUTOFF = float(sys.argv[4])
    print("# calculated by ts_calc_cv.py")
    print("# Search query: {:s}".format(sys.argv[1]))
    print("# nframe_skip: {:d}".format(nframe_skip))
    print("# nframe_end: {:d}".format(nframe_end))
    print("# E_CUTOFF: {:f}".format(E_CUTOFF))

    Ts = []
    e_sq = {}
    e = {}
    n = {}

    for tsfile in glob.glob(sys.argv[1]):
        ts = TsFile(tsfile)
        ts.open_to_read()
        ts.read_header()
        
        iframe = 0
        while ts.has_more_data():

            tsdata, tslines = ts.read_onestep()
            iframe += 1
            if iframe <= nframe_skip:
                continue
            if iframe > nframe_end:
                break

            T = float(tsdata[0][ts.head_col.temp])
            etot = float(tsdata[0][ts.head_col.e_tot])
            #evelo = float(tsdata[0][ts.head_col.e_velo])

            if etot > E_CUTOFF:
                continue

            if T in Ts:
                e_sq[T] += etot * etot
                e[T] += etot
                n[T] += 1
            else:
                Ts.append(T)
                e_sq[T] = etot * etot
                e[T] = etot
                n[T] = 1

    for T in sorted(Ts):
        fn = float(n[T])
        print('{:6.2f} {:8.3f} {:12d}'.format(T, ((e_sq[T]/fn - (e[T]/fn)**2) / (BOLTZC*T*T)), n[T]))

