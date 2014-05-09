#!/usr/bin/env python

import sys

if len(sys.argv) < 6 :
    print ('')
    print (' Usage: % SCRIPT [template tcl file] [eigen vector file] [(id begin, id end) ...] [output tcl file]')
    print ('')
    sys.exit(2)
    
# input for files
f_tcl_in = open(sys.argv[1], 'r')
f_pca_in = open(sys.argv[2], 'r')
f_tcl_out = open(sys.argv[-1], 'w')

# input for ID pairs
id_pairs = []
n = 1
for arg in sys.argv[3:-1] :
    if n == 1:
        tp = (int(arg),)
    else :
        id_pairs.append(tp + (int(arg),))
    n *= -1
if n == -1:
    print ('')
    print (' Usage: % SCRIPT [template tcl file] [eigen vector file] [(id begin, id end) ...] [output tcl file]')
    print ('')
    sys.exit(2)
    
ev = []
i_xyz = 0
for line in f_pca_in :
    if line.find('#') != -1 :
        continue
    i_xyz += 1
    if i_xyz == 1:
        vec = (float(line.strip()), )
    elif i_xyz == 2:
        vec += (float(line.strip()), )
    elif i_xyz == 3:
        vec += (float(line.strip()), )
        ev.append(vec)
        i_xyz = 0
        
for line_tcl in f_tcl_in :
    
    if line_tcl.find('##SELECTION##') != -1 :
        tcl = ''
        for (i,pair) in enumerate(id_pairs) :
            if i == 0:
                tcl += (' serial %i to %i' % pair)
            else :
                tcl += (' or serial %i to %i' % pair)
        f_tcl_out.write(line_tcl.replace('##SELECTION##', tcl))
        
    elif line_tcl.find('##VECTORS##') != -1 :
        tcl = ''
        for vec in ev :
            tcl += ' {%f %f %f}' % vec
        f_tcl_out.write(line_tcl.replace('##VECTORS##', tcl))
        
    else :
        f_tcl_out.write(line_tcl)
