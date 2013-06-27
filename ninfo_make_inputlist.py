#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2012/01/23
@author: Naoto Hori
'''

import sys
import glob
import re

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print 'Usage: SCRIPT [input DIR] [prefix] [output file]'
        sys.exit(2)
        
    files = glob.glob(sys.argv[1] + '/*.ninfo')
    prefix = sys.argv[2]
    max_unitnum = 0
    pairs = set()
    others = set()
    
    re_ddd_ddd = re.compile('.*/%s(\d\d\d)_(\d\d\d)\.ninfo' % prefix)
    re_others = re.compile('.*/([^/]*.ninfo)')
    for fname in files:
        reobj = re_ddd_ddd.match(fname)
        if reobj:
            unit1 = int(reobj.group(1))
            unit2 = int(reobj.group(2))
            if unit1 > max_unitnum:
                max_unitnum = unit1
            if unit2 > max_unitnum:
                max_unitnum = unit2
            pairs.add((unit1,unit2))
        else:
            reobj2 = re_others.match(fname)
            if reobj2:
                others.add(reobj2.group(1))
            else:
                print 'Error reobj2 is None'
                sys.exit(2)
            
    file_out = open(sys.argv[3],'w')
    n = 0
    for i in xrange(1,max_unitnum+1):
        for j in xrange(i,max_unitnum+1):
            if (i,j) in pairs:
                n += 1
                file_out.write('NINFO(%i/%i) %i\n' % (i,j,n))
                file_out.write('%i = %s%03i_%03i.ninfo\n' % (n,prefix,i,j))
    for other in others:
        file_out.write('NINFO(/)\n')
        file_out.write(' = ' + other+'\n')
        
        