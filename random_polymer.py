#!/usr/bin/env python

import numpy as np
import math
import NeRF
import random
import sys

def random_polymer(N, bond_length, bond_angle, Dexv):

    xyz = []
    
    # 1st particle at (0, 0, 0)
    xyz.append(np.array([0., 0., 0.]))
    
    # 2nd particle at (0, 0, L)
    xyz.append(np.array([0., 0., bond_length]))
    
    # 3rd particle at (Lsin(theta), 0, L + Lcos(theta))
    theta = math.pi - bond_angle
    xyz.append(np.array([bond_length*math.sin(theta), 0., bond_length*(1.0+math.cos(theta))]))
    
    i = 3
    while i < N:
    
        n_attempt = 0
        while True:
            p = random.uniform(0., 2*math.pi)
            new = NeRF.NeRF(xyz[-3], xyz[-2], xyz[-1], bond_length, bond_angle, p)
    
            flg_pass = True
            for j in range(0, i-1):
                if np.linalg.norm(new - xyz[j]) < Dexv:
                    flg_pass = False
                    break
        
            if flg_pass:
                xyz.append(new)
                i += 1
                break
    
            n_attempt += 1
            if n_attempt == 20:
                # If 20 attempts do not work out, go 3 particles back
                xyz.pop()
                xyz.pop()
                xyz.pop()
                i -= 3

    return xyz

if __name__ == "__main__":

    import sys

    if len(sys.argv) != 3:
        print('Usage: SCRIPT [input FASTA] [output XYZ]')
        sys.exit(0)

    #N = 5000
    L = 5.9  # bond length
    t = 2.618   # radian 
    Dexv = L

    # FASTA input
    seq = ''
    for l in open(sys.argv[1]):
        if l.startswith('>') or l.startswith('#'):
            continue
        if len(l.strip()) == 0:
            continue

        seq += l.strip().replace(' ','').replace('*','')
    
    N = len(seq)

    xyz = random_polymer(N, L, t, Dexv)

    #''' XYZ format '''
    f = open(sys.argv[2], 'w')
    for i in range(N):
        f.write('%s  %f %f %f\n' % (seq[i], xyz[i][0], xyz[i][1], xyz[i][2]))
    f.close()

