#!/usr/bin/env python
'''
@author: Naoto Hori
'''
import sys
import glob
import re

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print('Usage: %SCRIPT [input dir] [unit list file] [OUT input] [OUT list]')
        sys.exit(2)
        
    re_unitlist = re.compile('^(\d+)\s+(\S+)\s+(\S+)$')
    units = {}
    maxunit = 0
    for line in open(sys.argv[2],'r'):
        reobj = re_unitlist.match(line)
        if reobj:
            unitid = int(reobj.group(1))
            if unitid > maxunit:
                maxunit = unitid
            #unittype = re.group(2)
            unitname = reobj.group(3)
            units[unitid] = unitname

    se = re.compile("(\d+)_(\d+)\.ninfo$")
    filenames = {}
    files = glob.glob(sys.argv[1] + '*.ninfo')
    for filename in files:
        s = se.search(filename)
        i = int(s.group(1))
        j = int(s.group(2))
        filenames[(i,j)] = filename

    file_input = open(sys.argv[3],'w')
    n = 0
    for i in range(1,maxunit+1):
        for j in range(i,maxunit+1):
            if (i,j) in filenames:
                n += 1
                file_input.write('NINFO(%i/%i) %i\n' % (i,j,n))
                
    n = 0
    for i in range(1,maxunit+1):
        for j in range(i,maxunit+1):
            if (i,j) in filenames:
                n += 1
                file_input.write('%i %s\n' % (n, filenames[(i,j)]))
        
    file_list = open(sys.argv[4],'w')
    n = 0
    for i in range(1,maxunit+1):
        for j in range(i,maxunit+1):
            if (i,j) in filenames:
                n += 1
                file_list.write('%i %i %i %s %s\n' % (n,i,j,units[i],units[j]))



