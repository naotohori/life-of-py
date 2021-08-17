#!/usr/bin/env python

import sys

if __name__ == "__main__":

    import sys

    if len(sys.argv) != 3:
        print('Usage: SCRIPT [input FASTA] [output XYZ]')
        sys.exit(0)

    L = 5.9  # bond length

    # FASTA input
    seq = ''
    for l in open(sys.argv[1]):
        if l.startswith('>') or l.startswith('#'):
            continue
        if len(l.strip()) == 0:
            continue

        seq += l.strip().replace(' ','').replace('*','')
    
    N = len(seq)

    xyz = []
    
    for i in range(N):
        xyz.append([0., 0., i*L])

    #''' XYZ format '''
    f = open(sys.argv[2], 'w')
    f.write('%d\n' % N)
    f.write('\n')
    for i in range(N):
        f.write('%s  %f %f %f\n' % (seq[i], xyz[i][0], xyz[i][1], xyz[i][2]))
    f.close()

