#!/usr/bin/env python

import sys
from cafysis.file_io.ninfo import NinfoFile
from cafysis.elements.ninfo import NinfoSet

if len(sys.argv) != 3:
    print ('Usage: SCRIPT [input ninfo] [output ninfo]')
    sys.exit(2)

ninfo = NinfoFile(sys.argv[1])
ninfo.open_to_read()
ns = NinfoSet()
ninfo.read_all(ns)
ninfo.close()

ns.update_info()
ns.sort_by_mp()
ns.reassign_id()

ninfo = NinfoFile(sys.argv[-1])
ninfo.open_to_write()
ninfo.write_all(ns)
ninfo.close()
