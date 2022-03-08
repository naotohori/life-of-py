#!/usr/bin/env python


import sys
from lop.file_io.func_line2ninfo import line2basestackDT
from lop.file_io.func_ninfo2line import basestackDTdist2line
from lop.file_io.ninfo import NinfoFile
from lop.para.rnaAform import ARNA
from lop.para.rnaDT15 import DT15
from lop.para.rnaDT13 import DT13
from lop.para.rnaNHT19 import NHT19
from lop.elements.ninfo import NinfoSet, BondLength, BondAngle, BaseStackDT, HBondDT, TertiaryStackDT


if len(sys.argv) != 4:
    print('Usage: [old ninfo] [model] [new ninfo]')
    print('model has to be either of DT13, DT15, CHT18, NHT19')
    sys.exit(2)

f_out = open(sys.argv[-1],'w')

for l in open(sys.argv[1]):
    if not l.startswith('bs-dist'):
        f_out.write(l)
        continue
    
    info = line2basestackDT(l, fmt=1)
    ss = info.type[0] + info.type[2]

    if sys.argv[2] == 'DT13':
        info.h = DT13.ST_U0[ss][0]
        info.s = DT13.ST_U0[ss][1]
        info.Tm = DT13.ST_U0[ss][2]

    elif sys.argv[2] == 'DT15':
        info.h = DT15.ST_U0[ss][0]
        info.s = DT15.ST_U0[ss][1]
        info.Tm = DT15.ST_U0[ss][2]

    elif sys.argv[2] == 'CHT18':
        info.h = CHT18.ST_U0[ss][0]
        info.s = CHT18.ST_U0[ss][1]
        info.Tm = CHT18.ST_U0[ss][2]

    elif sys.argv[2] == 'NHT19':
        info.h = NHT19.ST_U0[ss][0]
        info.s = NHT19.ST_U0[ss][1]
        info.Tm = NHT19.ST_U0[ss][2]

    else:
        print('Error: unknown model')
        sys.exit(2)

    f_out.write( basestackDTdist2line(info) )

