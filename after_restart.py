#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2014/06/17
@author: Naoto Hori
'''

import sys
from cafysis.file_io.dcd import DcdFile
from cafysis.file_io.ts import TsFile
from copy import copy

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print 'Usage: % SCRIPT [name] [DCD1] [DCD2] ([DCD3] ...) [output dir]'
        sys.exit(2)
    
    dirs = sys.argv[2:-1]
    name = sys.argv[1]
    dir_out = sys.argv[-1]

    # Prepare output files
    out_ts = TsFiles( '%s/%s_%04i.ts' % (dir_out, prefix, i) )
    out_ts.open_to_write()

    out_dcd = DcdFiles( '%s/%s_%04i.dcd' % (dir_out, prefix, i) )
    out_dcd.open_to_write()

    first = True
    for d in dirs:
        in_ts = TsFiles( '%s/%s_%04i.ts' % (d, prefix, i) )
        in_ts.open_to_read()
        in_ts.read_header()
    
        in_dcd = DcdFiles( '%s/%s_%04i.dcd' % (d, prefix, i) )
        in_dcd.open_to_read()
        in_dcd.read_header()
    
        # For the first directory
        if first:
            out_ts.header_linse = in_ts.header_lines
            out_ts.write_header()
            out_dcd.set_header( in_dcd.get_header() )
            out_dcd.write_header()
            first = False
        else:
            in_ts.read_onestep()
            in_dcd.read_onestep()
    
        nstep_dcd = 0
        while in_dcd.has_more_data():
            out_dcd.write_onestep( in_dcd.read_onestep() )
            nstep_dcd += 1
    
        nstep_ts = 0
        while in_ts.has_more_data():
            tsdata, tslines = in_ts.read_onestep()
            out_ts.write_onestep( tslines )
            nstep_ts += 1
        
        if nstep_ts != nstep_dcd:
            print( 'nstep_ts != nstep_dcd' )
        
        in_ts.close()
        in_dcd.close()
    
