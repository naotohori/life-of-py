#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@author: Naoto Hori
'''

import sys
from lop.file_io.sisbp import SisbpFile

def sisbp_concatenate(filepaths):
    n_frames = []  # To return the number of frames

    # Number of input Sisbp files
    num_sisbp = len(filepaths) - 1

    # New Sisbp file
    f_out  = SisbpFile(filepaths[-1])
    f_out.open_to_write()

    f_in = SisbpFile(filepaths[0])
    f_in.open_to_read()
    f_in.read_header() 

    # Copy the header
    f_out.set_header(f_in.get_header())
    f_out.write_header()

    #print filepaths[0], f_in.get_header().nset
    #n_frames.append(f_in.get_header().nset)
    i = 0
    while f_in.has_more_data() :
        f_out.write_onestep(*f_in.read_onestep())
        i += 1
    n_frames.append(i)
    f_in.close()
    
    for i in range(1,num_sisbp) :
        f_in = SisbpFile(filepaths[i])
        f_in.open_to_read()
        f_in.read_header()
        #f_in.skip_onestep()  # skip the first step
        #print filepaths[i], f_in.get_header().nset - 1
        #n_frames.append(f_in.get_header().nset - 1)
        i = 0
        while f_in.has_more_data() :
            f_out.write_onestep(*f_in.read_onestep())
            i += 1
        n_frames.append(i)
        f_in.close()

    f_out.close()

    return n_frames

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print('Usage: % SCRIPT [DCD1] [DCD2] ([DCD3] ...) [output DCD]')
        sys.exit(2)
    n_frames = sisbp_concatenate(sys.argv[1:])
    for i,n in enumerate(n_frames):
        print('file ',i+1,': ', n)
    print('Total ', sum(n_frames))
