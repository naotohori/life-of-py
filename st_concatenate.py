#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2017/01/04
@author: Naoto Hori
'''

import sys

def st_concatenate(dir_in, dir_out, nrun):
   
    f_out = open('%s/azo.st' % dir_out,'w')
   
    ### Run 1
    f_in = open('%s/run001/azo.st' % dir_in,'r')
   
    # Write step=0
    f_out.write( f_in.readline() )

    # Skip step=1
    f_in.readline()

    # Read & write the rest of file
    line = f_in.readline()
    while(line):
        f_out.write( line )
        line = f_in.readline()
   
    f_in.close()
   
    ### Run 2 - nrun 
    for irun in range(2, nrun+1):
        f_in = open('%s/run%03i/azo.st' % (dir_in,irun,),'r')
   
        # Skip the first step (same as the last step of previous run)
        f_in.readline()

        # Read and write the rest of file
        line = f_in.readline()
        while(line):
            f_out.write( line )
            line = f_in.readline()
   
        f_in.close()
   
    f_out.close()

if __name__ == '__main__':
        
    if len(sys.argv) != 4:
        print 'Usage: % SCRIPT [input root dir] [output dir] [last run number]'
        sys.exit(2)

    dir_in = sys.argv[1]
    dir_out = sys.argv[2]
    nrun = int(sys.argv[3])

    if nrun < 3:
        print 'Usage: Error nrun < 3'
        sys.exit(2)

    st_concatenate(dir_in, dir_out, nrun)
