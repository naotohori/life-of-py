#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2014/06/17
@author: Naoto Hori
'''

import sys
from lop.file_io.dcd import DcdFile
from lop.file_io.ts import TsFile

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print('Usage: % SCRIPT [name] [Dir 1] ([Dir 2] [Dir 3] ...) [output dir]')
        sys.exit(2)
    
    name = sys.argv[1]
    dirs = sys.argv[2:-1]
    dir_out = sys.argv[-1]

    # Prepare output files
    out_ts = TsFile( '%s/%s.ts' % (dir_out, name) )
    out_ts.open_to_write()

    out_dcd = DcdFile( '%s/%s.dcd' % (dir_out, name) )
    out_dcd.open_to_write()

    nstep_total_ts = 0
    nstep_total_dcd = 0
    for (i_dir, d) in enumerate(dirs):
        in_ts = TsFile( '%s/%s.ts' % (d, name) )
        in_ts.open_to_read()
        in_ts.read_header()
    
        in_dcd = DcdFile( '%s/%s.dcd' % (d, name) )
        in_dcd.open_to_read()
        in_dcd.read_header()
    
        # For the first directory
        if i_dir == 0:
            out_ts.header_linse = in_ts.header_lines
            out_ts.copy_header( in_ts )
            out_ts.write_header()
            out_dcd.set_header( in_dcd.get_header() )
            out_dcd.write_header()
        else:
            # restartの場合は、最初のフレームを除外
            in_ts.read_onestep()
            in_dcd.read_onestep()
    
        nstep_dcd = 0
        while in_dcd.has_more_data():
            out_dcd.write_onestep( in_dcd.read_onestep() )
            nstep_dcd += 1
    
        nstep_ts = 0
        while in_ts.has_more_data():
            # step=1を除外(restartの時はない)
            if i_dir == 0 and nstep_ts == 1:
                in_ts.skip_onestep()
                nstep_ts += 1
                continue
            tsdata, tslines = in_ts.read_onestep()
            out_ts.write_onestep( tslines )
            nstep_ts += 1
        if i_dir == 0:
            # step=1を除外した分を引く
            nstep_ts -= 1
        
        if i_dir == 0:
            print(('TS  in %s: #frame = %i (step=1 is eliminated.)' % (d, nstep_ts)))
            print(('DCD in %s: #frame = %i' % (d, nstep_dcd)))
        else:
            print(('TS  in %s: #frame = %i (the first step is eliminated.)' % (d, nstep_ts)))
            print(('DCD in %s: #frame = %i (the first step is eliminated.)' % (d, nstep_dcd)))

        if nstep_ts != nstep_dcd:
            print( 'Warning: nstep_ts != nstep_dcd' )

        nstep_total_ts += nstep_ts
        nstep_total_dcd += nstep_dcd

        in_ts.close()
        in_dcd.close()
    
    print(('Total # of frames in  TS output to %s: %i' % (dir_out, nstep_total_ts)))
    print(('Total # of frames in DCD output to %s: %i' % (dir_out, nstep_total_dcd)))
