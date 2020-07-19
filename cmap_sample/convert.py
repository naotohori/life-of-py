#!/usr/bin/env python

import sys
import numpy as np

VAL_STRAND = 3
VAL_HELIX = 2
VAL_OTHER = 1
VAL_NO = 0


#if len(sys.argv) != 4:
#    print ('Usage: SCRIPT [input cor] [# contact (52)] [output gnudat]')
#    sys.exit(2)

#num_con = 196


####### Read dat file
if len(sys.argv) != 4:
    print ('Usage: SCRIPT [input dat (HeatMap)] [input dssp_struct] [output gnudat]')
    sys.exit(2)

file_path = sys.argv[1]
dssp_path = sys.argv[2]
outfile_path = sys.argv[-1]

num_res = 0
nline = 0
for l in open(file_path,'r'):
    if num_res == 0:
        num_res = len(l.split())
    nline += 1
if num_res != nline:
    print ('num_res != nline')
    sys.exit(2)


######## Read DSSP file
helix_res0 = []
strand_res0 = []
nline = 0
for l in open(dssp_path,'r'):
    lsp = l.split()
    i = int(lsp[0])
    nline += 1
    if len(lsp) != 2:
        continue
    c = lsp[1]
    if c in ('H','G','I'):
        helix_res0.append(i-1)
    elif c in ('E'):
        strand_res0.append(i-1)
if num_res != nline:
    print ('nline(dssp) != num_res')
    sys.exit(2)

    
dat = np.zeros((num_res+1, num_res+1))
for i,l in enumerate(open(file_path,'r')):
    #if l.strip() == '':
    #    continue
    lsp = l.split()
    for j in range(num_res):
        if int(lsp[j]) == 1:
            if i in helix_res0 and j in helix_res0:
                dat[i+1,j+1] = VAL_HELIX
                dat[j+1,i+1] = VAL_HELIX
            elif i in strand_res0 and j in strand_res0:
                dat[i+1,j+1] = VAL_STRAND
                dat[j+1,i+1] = VAL_STRAND
            else:
                dat[i+1,j+1] = VAL_OTHER
                dat[j+1,i+1] = VAL_OTHER
        else:
            dat[i+1,j+1] = VAL_NO
            dat[j+1,i+1] = VAL_NO


f_out = open(outfile_path,'w')

f_out.write('%i %i %f\n' % (0,0,dat[1,1]))
for j in range(1,num_res):
    f_out.write('%i %i %f\n' % (0,j,dat[1,j]))
    f_out.write('%i %i %f\n' % (0,j,dat[1,j+1]))
f_out.write('%i %i %f\n' % (0,num_res,dat[1,num_res]))
f_out.write('\n')

for i in range(1,num_res):
    f_out.write('%i %i %f\n' % (i,0,dat[i,1]))
    for j in range(1,num_res):
        f_out.write('%i %i %f\n' % (i,j,dat[i,j]))
        f_out.write('%i %i %f\n' % (i,j,dat[i,j+1]))
    f_out.write('%i %i %f\n' % (i,num_res,dat[i,num_res]))
    f_out.write('\n')

    f_out.write('%i %i %f\n' % (i,0,dat[i+1,1]))
    for j in range(1,num_res):
        f_out.write('%i %i %f\n' % (i,j,dat[i+1,j]))
        f_out.write('%i %i %f\n' % (i,j,dat[i+1,j+1]))
    f_out.write('%i %i %f\n' % (i,num_res,dat[i+1,num_res]))
    f_out.write('\n')

f_out.write('%i %i %f\n' % (num_res,0,dat[num_res,1]))
for j in range(1,num_res):
    f_out.write('%i %i %f\n' % (num_res,j,dat[num_res,j]))
    f_out.write('%i %i %f\n' % (num_res,j,dat[num_res,j+1]))
f_out.write('%i %i %f\n' % (num_res,num_res,dat[num_res,num_res]))
f_out.write('\n')
