#!/usr/bin/env python
#vim:fileencoding=UTF-8

from cafysis.file_io.dcd import DcdFile
import copy
import sys

if len(sys.argv) != 3:
    print 'Usage: SCRIPT [input DCD] [output DCD]'
    sys.exit(2)

#ID_DOM_INI = int(sys.argv[2]) - 1  # 重心を求める際に必要
#ID_DOM_END = int(sys.argv[3]) - 1

dcd = DcdFile(sys.argv[1])
dcd.open_to_read()
dcd_out = DcdFile(sys.argv[-1])
dcd_out.open_to_write()

# header
dcd.read_header()
header = dcd.get_header()
header_new = copy.deepcopy(header)
header_new.nunit_real = 1
header_new.nmp_real = 196
header_new.lunit2mp = []
header_new.lunit2mp.append(196)
dcd_out.set_header(header_new)
dcd_out.write_header()

nmp = header.nmp_real

while dcd.has_more_data() :

    data = dcd.read_onestep()

    data_S = []
    for i in range(0,587,3):
        data_S.append(data[i])
    
    dcd_out.write_onestep(data_S) 
    
dcd.close()
dcd_out.close()
