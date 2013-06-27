#!/usr/bin/env python

file_out = open('edc_pre.3i8h.out','w')
import sys

if len(sys.argv) != 3:
    print ('Usage: SCRIPT [input (edc_pre.residue)] [output (edc_pre.residue.mpunit)]')
    sys.exit(2)
    
file_out = open(sys.argv[2], 'w')

ichain_now = 0
ires_offset = 0
imp_offset = 0
ires_pre = 0
imp_pre  = 0
for line in open(sys.argv[1],'r'):
    if line.find('#') != -1:
        file_out.write(line)
        continue
    linesp = line.split()
    ichain = int(linesp[0])
    ires   = int(linesp[1])
    imp    = int(linesp[2])
    ires_l = int(linesp[3])
    imp_l  = int(linesp[4])
    if ichain != ichain_now:
        ichain_now = ichain
        ires_offset = ires_pre
        imp_offset  = imp_pre
    file_out.write(line.rstrip()+('  %5i %5i\n' % (ires - ires_offset, imp - imp_offset)))
    ires_pre = ires
    imp_pre  = imp

