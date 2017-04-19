#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2017/01/04
@author: Naoto Hori
'''

import sys

def ts_concatenate(dir_in, dir_out, nrun):
   
    f_out = open('%s/vpk.ts' % dir_out,'w')
   
    ## Write Header
    #f_out.write('#                step    tempk     radg       etot      velet qscore     rmsd\n')
    f_out.write('#                step    tempk     radg       etot      velet qscore     rmsd      local         go      repul      stack     hbond     elect\n')
   
   
    ### Run 1
    f_in = open('%s/run001/vpk.ts' % dir_in,'r')
    # Skip 9
    for i in range(9):
        f_in.readline()
   
    # Step 0
    f_out.write( f_in.readline() )
     
    # Skip step 1
    for i in range(1):
        f_in.readline()
   
    # Read and write the rest of file
    line = f_in.readline()
    while(line):
        f_out.write( line )
        line = f_in.readline()
   
    f_in.close()
   
    ### Run 2 - nrun 
    for irun in range(2, nrun+1):
        f_in = open('%s/run%03i/vpk.ts' % (dir_in,irun,),'r')
     
        # Skip 9
        for i in range(9):
            f_in.readline()
     
        # Skip the first step (same as the last step of previous run)
        for i in range(1):
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

    ts_concatenate(dir_in, dir_out, nrun)
