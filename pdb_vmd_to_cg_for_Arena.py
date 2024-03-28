#!/usr/bin/env python

def convert(filepath, outfilepath):

    fout = open(outfilepath, 'w')
    for l in open(filepath):
        if l.startswith('ATOM  '):
            #nuc = l[17:18]
            ##print (nuc)
            #newl = l[0:17] + '  ' + l[17:18] + l[20:] 

            nuc = l[13:14]
            newl = l[0:21] + nuc + l[22:] 

            fout.write(newl)
        else:
            fout.write(l)

if __name__ == "__main__":
    import glob
    import os.path

    import sys

    if len(sys.argv) == 1:
        print ('Usage: SCRIPT (dir)')
        print ('  or   SCRIPT (pdb) [(pdb) [(pdb)]...]')
        sys.exit(0)

    elif len(sys.argv) == 2:
        if os.path.isdir(sys.argv[1]):
            pdbfiles = glob.glob('./*.pdb')
        else:
            pdbfiles = [sys.argv[1], ]

    else:
        pdbfiles = sys.argv[1:]

    for pdb in pdbfiles:
        bn = os.path.basename(pdb)
        out = bn[0:-4] + '.cg.pdb'

        convert(pdb, out)
