#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2011/11/09
@author: Naoto Hori
'''

import sys
from cafysis.file_io.dcd import DcdFile
from copy import copy

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print('Usage: % SCRIPT [DCD1] [DCD2] ([DCD3] ...) [output DCD]')
        sys.exit(2)
        
    # Number of input DCD files
    num_dcd = len(sys.argv) - 2
    
    # New DCD file
    f_out  = DcdFile(sys.argv[-1])
    f_out.open_to_write()
    
    # Count the total frame number
    num_frame = 1
    for i in range(1,num_dcd) :
        f_in = DcdFile(sys.argv[i])
        f_in.open_to_read()
        f_in.read_header()
        num_frame += f_in.get_header().nset
        f_in.close()
        
    # Get the total step number from final DCD file
    f_in = DcdFile(sys.argv[num_dcd])
    f_in.open_to_read()
    f_in.read_header()
    num_frame += f_in.get_header().nset
    num_step = f_in.get_header().nstep
    f_in.close()
        
    f_in = DcdFile(sys.argv[1])
    f_in.open_to_read()
    f_in.read_header() 
    header = copy(f_in.get_header())
    header.nset = num_frame
    header.nstep = num_step
    f_out.set_header(header)
    f_out.write_header()
    print(sys.argv[1], f_in.get_header().nset)
    while f_in.has_more_data() :
        f_out.write_onestep(f_in.read_onestep())
    f_in.close()
    
    for i in range(2,num_dcd+1) :
        f_in = DcdFile(sys.argv[i])
        f_in.open_to_read()
        f_in.read_header()
        #f_in.skip_onestep()  # skip the first step
        print(sys.argv[i], f_in.get_header().nset)
        while f_in.has_more_data() :
            f_out.write_onestep(f_in.read_onestep())
        f_in.close()

    f_out.close()
