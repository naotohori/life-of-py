#!/usr/bin/env python

def convert(filepath, outfilepath):

    fout = open(outfilepath, 'w')
    for l in open(filepath):
        if l.startswith('ATOM  '):
            nuc = l[17:18]
            print (nuc)
            newl = l[0:17] + '  ' + l[17:18] + l[20:] 
            fout.write(newl)
        else:
            fout.write(l)

if __name__ == "__main__":
    import glob
    import os.path

    pdbfiles = glob.glob('./*.pdb')
    for pdb in pdbfiles:
        bn = os.path.basename(pdb)
        out = bn[0:-4] + '.cg.pdb'

        convert(pdb, out)
