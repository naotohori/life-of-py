#!/usr/bin/env python
#vim:fileencoding=UTF-8

from cafysis.file_io.dcd import DcdFile
import copy
import sys

if len(sys.argv) != 4:
    print 'Usage: SCRIPT [input DCD] [last index of RNA (count from 0)] [output DCD]'
    sys.exit(2)

last_mp = int(sys.argv[2])

dcd = DcdFile(sys.argv[1])
dcd.open_to_read()
dcd_out = DcdFile(sys.argv[-1])
dcd_out.open_to_write()

# header
nres = (last_mp - last_mp % 3) / 3 + 1
dcd.read_header()
header = dcd.get_header()
header_new = copy.deepcopy(header)
header_new.nunit_real = 1
header_new.nmp_real = nres
header_new.lunit2mp = []
header_new.lunit2mp.append(nres)
dcd_out.set_header(header_new)
dcd_out.write_header()

nmp = header.nmp_real

while dcd.has_more_data() :

    data = dcd.read_onestep()

    data_S = []
    for i in range(0,last_mp+1,3):
        data_S.append(data[i])
    
    dcd_out.write_onestep(data_S) 
    
dcd.close()
dcd_out.close()
