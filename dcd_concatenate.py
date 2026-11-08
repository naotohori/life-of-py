#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2011/11/09
@author: Naoto Hori
'''

import argparse
from lop.file_io.dcd import DcdFile
from lop.dcd_frame_count import count_frame
from copy import copy
        

def dcd_concatenate(filepaths, skip_first=True):
    '''Concatenate DCD files. filepaths = [input1, input2, ..., output].
    If skip_first is True, the first frame of the second and later inputs
    is skipped (it duplicates the last frame of the previous run on restart).'''
    n_frames = []  # To return the number of frames

    # Number of input DCD files
    num_dcd = len(filepaths) - 1
    
    # New DCD file
    f_out  = DcdFile(filepaths[-1])
    f_out.open_to_write()
    
    # Count the total frame number
    num_frame = 1 if skip_first else 0
    #for i in range(0,num_dcd-1) :
    for i in range(0,num_dcd) :
        #f_in = DcdFile(filepaths[i])
        #f_in.open_to_read()
        #f_in.read_header()
        #nset = f_in.get_header().nset
        #if nset is not None:
        #    num_frame += nset - 1
        #f_in.close()
        if skip_first:
            num_frame += count_frame(filepaths[i]) - 1
        else:
            num_frame += count_frame(filepaths[i])
        
    # Get the total step number from final DCD file
    f_in = DcdFile(filepaths[num_dcd-1])
    f_in.open_to_read()
    f_in.read_header()
    #nset = f_in.get_header().nset
    #if nset is not None:
    #    num_frame += nset - 1
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
    #n_frames.append(f_in.get_header().nset)
    i = 0
    while f_in.has_more_data() :
        coords = f_in.read_onestep()
        if f_in.get_header().with_unit_cell:
            f_out.get_header().unit_cell_xyz = f_in.get_header().unit_cell_xyz
            f_out.get_header().unit_cell_abc = f_in.get_header().unit_cell_abc
        f_out.write_onestep(coords)
        i += 1
    n_frames.append(i)
    f_in.close()
    
    for i in range(1,num_dcd) :
        f_in = DcdFile(filepaths[i])
        f_in.open_to_read()
        f_in.read_header()
        if skip_first:
            f_in.skip_onestep()  # skip the first step
        #print filepaths[i], f_in.get_header().nset - 1
        #n_frames.append(f_in.get_header().nset - 1)
        i = 0
        while f_in.has_more_data() :
            coords = f_in.read_onestep()
            if f_in.get_header().with_unit_cell:
                f_out.get_header().unit_cell_xyz = f_in.get_header().unit_cell_xyz
                f_out.get_header().unit_cell_abc = f_in.get_header().unit_cell_abc
            f_out.write_onestep(coords)
            i += 1
        n_frames.append(i)
        f_in.close()

    f_out.close()

    return n_frames


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Concatenate DCD files. By default, the first frame of '
                    'the second and later inputs is skipped, as it duplicates '
                    'the last frame of the previous run (restart).')
    parser.add_argument('inputs', nargs='+', help='input DCD files (two or more)')
    parser.add_argument('output', help='output DCD file')
    parser.add_argument('--noskip', action='store_true',
                        help='keep the first frame of every input DCD')
    args = parser.parse_args()
    if len(args.inputs) < 2:
        parser.error('at least two input DCD files are required')
    n_frames = dcd_concatenate(args.inputs + [args.output],
                               skip_first=not args.noskip)
    for i,n in enumerate(n_frames):
        print('file ',i+1,': ', n)
