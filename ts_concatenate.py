#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2017/01/04
@author: Naoto Hori
'''

import sys

def ts_concatenate(dir_in, dir_out, nrun):
   
    f_out = open('%s/azo.ts' % dir_out,'w')
    f_out_rna = open('%s/azo.rna.ts' % dir_out,'w')
    f_out_ion = open('%s/azo.ion.ts' % dir_out,'w')
    f_out_rnaion = open('%s/azo.rna-ion.ts' % dir_out,'w')
   
    ## Write Header
    #f_out.write('#                step    tempk     radg       etot      velet qscore     rmsd\n')
    f_out.write('#                step    tempk     radg       etot      velet qscore     rmsd      local         go      repul      stack     tstack      hbond     thbond      elect\n')
    f_out_rna.write('#                step    tempk     radg       etot      velet qscore     rmsd      local         go      repul      stack     tstack      hbond     thbond      elect\n')
    f_out_ion.write('#                step    tempk     radg       etot      velet qscore     rmsd      local         go      repul      stack     tstack      hbond     thbond      elect\n')
    f_out_rnaion.write('#                step    tempk     radg       etot      velet qscore     rmsd      local         go      repul      stack     tstack      hbond     thbond      elect\n')
   
   
    ### Run 1
    f_in = open('%s/run001/azo.ts' % dir_in,'r')
    # Skip 9
    for i in range(9):
        f_in.readline()
   
    # Step 0
    f_out.write( f_in.readline() )
    f_out_rna.write( f_in.readline().replace('#1','  ') )
    f_out_ion.write( f_in.readline().replace('#2','  ') )
    f_out_rnaion.write( f_in.readline().replace('#1   -2','       ') )
     
    # Skip step 1
    for i in range(4):
        f_in.readline()
   
    # Read and write the rest of file
    line = f_in.readline()
    iline = 0
    while(line):
        if iline == 0:
            f_out.write( line )
        elif iline == 1:
            f_out_rna.write( line.replace('#1','  ') )
        elif iline == 2:
            f_out_ion.write( line.replace('#2','  ') )
        elif iline == 3:
            f_out_rnaion.write( line.replace('#1   -2','       ') )
   
        line = f_in.readline()
        iline += 1
        if iline == 4:
            iline = 0
   
    if iline != 0:
        print 'Error: iline != 0, irun = ',irun
     
    f_in.close()
   
   
    ### Run 2 - nrun 
    for irun in range(2, nrun+1):
        f_in = open('%s/run%03i/azo.ts' % (dir_in,irun,),'r')
     
        # Skip 9
        for i in range(9):
            f_in.readline()
     
        # Skip the first step (same as the last step of previous run)
        for i in range(4):
            f_in.readline()
     
        # Read and write the rest of file
        line = f_in.readline()
        iline = 0
        while(line):
            if iline == 0:
                f_out.write( line )
            elif iline == 1:
                f_out_rna.write( line.replace('#1','  ') )
            elif iline == 2:
                f_out_ion.write( line.replace('#2','  ') )
            elif iline == 3:
                f_out_rnaion.write( line.replace('#1   -2','       ') )
     
            line = f_in.readline()
            iline += 1
            if iline == 4:
                iline = 0
     
        if iline != 0:
            print 'Error: iline != 0, irun = ',irun
     
        f_in.close()
     
    f_out.close()
    f_out_rna.close()
    f_out_ion.close()
    f_out_rnaion.close()
    
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
