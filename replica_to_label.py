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
    if len(sys.argv) < 5:
        print('Usage: %SCRIPT [input Dir] [name] [final ID] [output DIR]')
        print('       ID shoule be from 0 to final ID.')
        sys.exit(2)

    dir_in = sys.argv[1]
    name = sys.argv[2]
    id_end = int(sys.argv[3])
    dir_out = sys.argv[-1]

    out_ts_files = []
    out_dcd_files = []
    #Prepare output files
    for id_rep in range(1, id_end+1):
        ts = TsFile('%s/%s_%04i.ts' % (dir_out, name, id_rep))
        ts.open_to_write()
        out_ts_files.append(ts)

        dcd = DcdFile('%s/%s_%04i.dcd' % (dir_out, name, id_rep))
        dcd.open_to_write()
        out_dcd_files.append(dcd)

    in_ts_files = []
    in_dcd_files = []
    #Open input files
    for id_rep in range(1, id_end+1):
        idx = id_rep - 1

        ts = TsFile('%s/%s_%04i.ts' % (dir_in, name, id_rep))
        ts.open_to_read()
        ts.read_header()
        in_ts_files.append(ts)

        out_ts_files[idx].copy_header(ts)
        out_ts_files[idx].write_header()

        dcd = DcdFile('%s/%s_%04i.dcd' % (dir_in, name, id_rep))
        dcd.open_to_read()
        dcd.read_header()
        in_dcd_files.append(dcd)

        out_dcd_files[idx].set_header( dcd.get_header() )
        out_dcd_files[idx].write_header()

    while in_dcd_files[0].has_more_data():
        for idx in range(id_end):
            (tsdata, tslines) = in_ts_files[idx].read_onestep()
            coord_matrix = in_dcd_files[idx].read_onestep()

            step  = int( tsdata[0][in_ts_files[idx].head_col.step] )
            label = int( tsdata[0][in_ts_files[idx].head_col.label] )
            idx_out = label - 1

            if step == 1:
                (tsdata, tslines) = in_ts_files[idx].read_onestep()
                step  = int( tsdata[0][in_ts_files[idx].head_col.step] )
                label = int( tsdata[0][in_ts_files[idx].head_col.label] )
                idx_out = label - 1

            out_ts_files[idx_out].write_onestep( tslines )
            out_dcd_files[idx_out].write_onestep( coord_matrix )


'''
    # For the last one step
    for idx in range(id_end):
        (tsdata, tslines) = in_ts_files[idx].read_onestep()
        coord_matrix = in_dcd_files[idx].read_onestep()

        label = int( tsdata[0][in_ts_files[idx].head_col.label] )
        idx_out = label - 1

        out_ts_files[idx_out].write_onestep( tslines )
        out_dcd_files[idx_out].write_onestep( coord_matrix )
'''
