#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2012/01/19
@author: Naoto Hori
'''

import sys
from cafysis.file_io.ninfo import NinfoFile
from cafysis.elements.ninfo import NinfoSet

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print 'Usage: %SCRIPT [input ninfo file] [output prefix]'
        sys.exit(2)
    
    file_ninfo = NinfoFile(sys.argv[1])
    file_ninfo.open_to_read()
    ns = NinfoSet()
    file_ninfo.read_all(ns)
    ns.update_info()
    file_ninfo.close()
    
    prefix = sys.argv[2]
    
    ndi = ns.dict_of_ninfoset_by_unit()
    for i in xrange(1,ns.max_unit+1):
        for j in xrange(i,ns.max_unit+1):
            if ndi[(i,j)].max_unit != 0:
                file_out = NinfoFile(prefix+'%03i_%03i.ninfo'%(i,j))
                file_out.open_to_write()
                file_out.write_all(ndi[(i,j)])
                #file_out.write_unit(ns, i, j)
                file_out.close()
    