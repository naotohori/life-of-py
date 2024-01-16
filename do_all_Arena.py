#!/usr/bin/env python

if __name__ == "__main__":
    import glob
    import os.path
    import subprocess

    pdbfiles = glob.glob('./*.cg.pdb')
    for pdb in pdbfiles:
        bn = os.path.basename(pdb)
        out = bn[0:-7] + '.aa.pdb'

        subprocess.call(['Arena', pdb, out])
