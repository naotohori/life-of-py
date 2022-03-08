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

from lop.file_io.pdb import PdbFile
from lop.elements.coord import Coord

if len(sys.argv) not in (2,3):
    print('Usage: SCRIPT [PDB]  (score cutoff = 0.5)')
    print('  or : SCRIPT [PDB] [score cutoff]')
    sys.exit(2)

SCORE_CUT = 0.5
if len(sys.argv) > 2:
    SCORE_CUT = float(sys.argv[2])

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
            print('Error: can not find atom "a" in ', ic+1, '-th chain, ',ir+1, '-th residue, ', ntd)
            sys.exit(2)
        if xyz_b is None:
            print('Error: can not find atom "b" in ', ic+1, '-th chain, ',ir+1, '-th residue, ', ntd)
            sys.exit(2)
        
        xyz_a -= com
        xyz_b -= com

        normal = xyz_a.cross(xyz_b) 
        normals.append(normal)

        locations.append( (ic+1,ir+1) )

factor_r3 = 1.0 / ((5.0-3.5)/(3.5**3))   # To normalize the function for distance (~28.58333)

def stack_score(comi, comj, normali, normalj):

        normali_abs = normali.norm()
        normalj_abs = normalj.norm()

        score = 0.0

        # Distance
        d = comi.distance(comj)
        if d > 5.0:
            #pass
            #continue
            return None, None, None
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
            #continue
            return None, None, None
        elif omega <= 25.0:
            score += 1.0
        else:
            score += 2.0 - 0.04 * omega

        #if score < 0.5:
        #    continue

        # Xi
        xi = math.asin( normali.cross(normalj).norm() / (normali.norm() * normalj.norm()) ) * 180.0 / math.pi
        #xi2 = math.asin( normali.cross(normalj*(-1)).norm() / (normali.norm() * normalj.norm()) ) * 180.0 / math.pi

        if normali.dot( v_ij ) > 0.0:
            sgni = -1
        else:
            sgni = +1
        if normalj.dot( v_ij ) > 0.0:
            sgnj = +1
        else:
            sgnj = -1
        
        if 45.0 < xi and xi < 135.0:
            #continue
            return None, None, None

        return score, (sgni, sgnj), (d, omega, xi)

''' Secondary stacking '''
for i in range(nres-1):
    comi = coms[i]
    normali = normals[i]

    comj = coms[i+1]
    normalj = normals[i+1]

    score, sgn, param = stack_score(comi, comj, normali, normalj)
    if score is not None and score >= SCORE_CUT:
        si, sj = sgn
        d, omega, xi = param
        print(('%s %+4i %+4i  %5.3f  %4.2f %6.2f %6.2f' % 
               ('S', si*locations[i][1], sj*locations[i+1][1], score, d, omega, xi)))

''' Tertiary stacking '''
for i in range(nres):
    comi = coms[i]
    normali = normals[i]

    for j in range(i+2, nres):
        comj = coms[j]
        normalj = normals[j]

        score, sgn, param = stack_score(comi, comj, normali, normalj)
        if score is not None and score >= SCORE_CUT:
            si, sj = sgn
            d, omega, xi = param
            print(('%s %+4i %+4i  %5.2f  %4.2f %6.2f %6.2f' % 
                   ('T', si*locations[i][1], sj*locations[j][1], score, d, omega, xi)))
