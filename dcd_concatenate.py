#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2011/11/09
@author: Naoto Hori
'''

import sys
from lop.file_io.dcd import DcdFile
from copy import copy
        

def dcd_concatenate(filepaths):
    n_frames = []  # To return the number of frames

    # Number of input DCD files
    num_dcd = len(filepaths) - 1
    
    # New DCD file
    f_out  = DcdFile(filepaths[-1])
    f_out.open_to_write()
    
    # Count the total frame number
    num_frame = 1
    for i in range(0,num_dcd-1) :
        f_in = DcdFile(filepaths[i])
        f_in.open_to_read()
        f_in.read_header()
        num_frame += f_in.get_header().nset - 1
        f_in.close()
        
    # Get the total step number from final DCD file
    f_in = DcdFile(filepaths[num_dcd-1])
    f_in.open_to_read()
    f_in.read_header()
    num_frame += f_in.get_header().nset - 1
    num_step = f_in.get_header().nstep
    f_in.close()
        
    f_in = DcdFile(filepaths[0])
    f_in.open_to_read()
    f_in.read_header() 
    header = copy(f_in.get_header())
    header.nset = num_frame
    header.nstep = num_step
    f_out.set_header(header)
    f_out.write_header()
    #print filepaths[0], f_in.get_header().nset
    n_frames.append(f_in.get_header().nset)
    while f_in.has_more_data() :
        f_out.write_onestep(f_in.read_onestep())
    f_in.close()
    
    for i in range(1,num_dcd) :
        f_in = DcdFile(filepaths[i])
        f_in.open_to_read()
        f_in.read_header()
        f_in.skip_onestep()  # skip the first step
        #print filepaths[i], f_in.get_header().nset - 1
        n_frames.append(f_in.get_header().nset - 1)
        while f_in.has_more_data() :
            f_out.write_onestep(f_in.read_onestep())
        f_in.close()

    f_out.close()

    return n_frames


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print('Usage: % SCRIPT [DCD1] [DCD2] ([DCD3] ...) [output DCD]')
        sys.exit(2)
    n_frames = dcd_concatenate(sys.argv[1:])
    for i,n in enumerate(n_frames):
        print('file ',i+1,': ', n)
