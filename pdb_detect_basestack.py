#!/usr/bin/env python
'''
Created on 2016/09/30
@author: Naoto Hori

### Ref. D.E. Condon and D.H. Turner et al. 
###      J. Chem. Theory Comput. (2015) 11: 2729
###      doi:10.1021/ct501025q

### Special thanks to Hung T. Nguyen
'''

import sys
import math

from cafysis.file_io.pdb import PdbFile
from cafysis.elements.coord import Coord

pf = PdbFile(sys.argv[1])
pf.open_to_read()
chains = pf.read_all()
pf.close()

coms = []
normals = []
locations = []

## Suppose all residues have a base

nres = 0
for ic, c in enumerate(chains):
    for ir, r in enumerate(c.residues):

        nres += 1

        ## A,U,G or C
        ntd = r.atoms[0].res_name.strip()

        ## Center of Mass, and find atoms "a" and "b"
        com = Coord()
        natom = 0
        xyz_a = None
        xyz_b = None
        for a in r.atoms:
            name = a.name.strip()
            if name.find("'") != -1:   # Skip ribose
                continue
            if name[0] == 'H': # Skip hydrogens
                continue
            if name.find("P") != -1:   # Skip phosphate
                continue
            natom += 1
            com += a.xyz

            if ntd == 'A':
                if name == 'N6':
                    xyz_a = a.xyz
                if name == 'C8':
                    xyz_b = a.xyz
            elif ntd == 'G':
                if name == 'O6':
                    xyz_a = a.xyz
                if name == 'C8':
                    xyz_b = a.xyz
            elif ntd == 'C':
                if name == 'O2':
                    xyz_a = a.xyz
                if name == 'N4':
                    xyz_b = a.xyz
            elif ntd == 'U':
                if name == 'O2':
                    xyz_a = a.xyz
                if name == 'O4':
                    xyz_b = a.xyz
    
        com /= float(natom)
        coms.append(com)

        if xyz_a is None:
            print 'Error: can not find atom "a" in ', ic+1, '-th chain, ',ir+1, '-th residue, ', ntd
            sys.exit(2)
        if xyz_b is None:
            print 'Error: can not find atom "b" in ', ic+1, '-th chain, ',ir+1, '-th residue, ', ntd
            sys.exit(2)
        
        xyz_a -= com
        xyz_b -= com

        normal = xyz_a.cross(xyz_b) 
        normals.append(normal)

        locations.append( (ic+1,ir+1) )

factor_r3 = 1.0 / ((5.0-3.5)/(3.5**3))   # To normalize the function for distance (~28.58333)

for i in range(nres):
    comi = coms[i]
    normali = normals[i]
    normali_abs = normali.norm()

    for j in range(i+2, nres):
        comj = coms[j]
        normalj = normals[j]
        normalj_abs = normalj.norm()

        score = 0.0

        # Distance
        d = comi.distance(comj)
        if d > 5.0:
            #print 11+locations[i][1], 11+locations[j][1], 'd=',d
            #continue
            pass
        elif d <= 3.5:
            score += 1.0
        else:
            score += factor_r3 * (5.0-d)/(d**3)

        # Omega
        v_ij = comi - comj
        omega = math.acos( normali.dot( v_ij ) / ( normali_abs * v_ij.norm() ) ) * 180.0 / math.pi
        if omega > 90.0:
            omega = 180.0 - omega

        if omega > 50.0:
            print 11+locations[i][1], 11+locations[j][1], 'omega=',omega
            continue
        elif omega <= 25.0:
            score += 1.0
        else:
            score += 2.0 - 0.04 * omega

        #if score < 0.5:
        #    continue

        # Xi
        xi = math.asin( normali.cross(normalj).norm() / (normali.norm() * normalj.norm()) ) * 180.0 / math.pi
        #xi2 = math.asin( normali.cross(normalj*(-1)).norm() / (normali.norm() * normalj.norm()) ) * 180.0 / math.pi
        print 11+locations[i][1], 11+locations[j][1], xi, score

