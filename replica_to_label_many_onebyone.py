#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2014/06/17
@author: Naoto Hori
'''

import sys
from cafysis.file_io.dcd import DcdFile
from cafysis.file_io.ts import TsFile
#import resource

#resource.setrlimit(resource.RLIMIT_NOFILE, (2000,2000))

if __name__ == '__main__':
    if len(sys.argv) < 5:
        print 'Usage: %SCRIPT [input Dir] [name] [final ID] [output DIR]'
        print '       ID shoule be from 0 to final ID.'
        sys.exit(2)

    dir_in = sys.argv[1]
    name = sys.argv[2]
    id_end = int(sys.argv[3])
    dir_out = sys.argv[-1]

    #Open input files
    in_ts_files = []
    in_dcd_files = []
    for id_rep in range(1, id_end+1):
        idx = id_rep - 1

        ts = TsFile('%s/%s_%04i.ts' % (dir_in, name, id_rep))
        ts.open_to_read()
        ts.read_header()
        in_ts_files.append(ts)

        dcd = DcdFile('%s/%s_%04i.dcd' % (dir_in, name, id_rep))
        dcd.open_to_read()
        dcd.read_header()
        in_dcd_files.append(dcd)

    # Loop for output label
    for id_lab in range(1, id_end+1):
        #Prepare output files
        tsout = TsFile('%s/%s_%04i.ts' % (dir_out, name, id_lab))
        tsout.open_to_write()

        dcdout = DcdFile('%s/%s_%04i.dcd' % (dir_out, name, id_lab))
        dcdout.open_to_write()

        #Rewind
        for id_rep in range(1, id_end+1):
            idx = id_rep - 1
            in_ts_files[idx].read_header()
            in_dcd_files[idx].read_header()

            if (id_rep == id_lab):
                tsout.copy_header(in_ts_files[idx])
                tsout.write_header()
                dcdout.set_header( in_dcd_files[idx].get_header() )
                dcdout.write_header()

        while in_dcd_files[0].has_more_data():
            for idx in range(id_end):
                (tsdata, tslines) = in_ts_files[idx].read_onestep()
                coord_matrix = in_dcd_files[idx].read_onestep()
    
                label = int( tsdata[0][in_ts_files[idx].head_col.label] )

                if (label == id_lab):
                    tsout.write_onestep( tslines )
                    dcdout.write_onestep( coord_matrix )
