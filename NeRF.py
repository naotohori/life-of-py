#!/usr/bin/env python

import numpy as np
from torsion import torsion

def NeRF(A, B, C, bond, angl, dihd):

    D2 = np.array([bond * np.cos(angl), 
                   bond * np.cos(dihd) * np.sin(angl), 
                   bond * np.sin(dihd) * np.sin(angl)])
  
    AB = B - A
    BC = C - B
    n_bc = BC / np.linalg.norm(BC)

    n = np.cross(AB, n_bc)
    n = n / np.linalg.norm(n)

    Mx = n_bc
    My = np.cross(n, n_bc)
    #My = My / np.linalg.norm(My)   # Not necessary as the norm should be 1
    Mz = n
    M = np.array([Mx, My, Mz]).T
    # Apply transpose becasue vectors (Mx, My, Mz) have to be the columns of M.

    D = np.matmul(M, D2) + C

    return D

def bond_length(A, B):
    AB = B - A
    return np.linalg.norm(AB)

def bond_angle(A, B, C):
    BA = A - B
    BC = C - B
    return np.arccos( np.dot(BA, BC) / (np.linalg.norm(BA) * np.linalg.norm(BC)) )


""" 
A sample data: coordinates of four particles
(in XYZ format (copy & paste to see in VMD) including the number of atoms and a blank line)
4

P      21.436  33.390  24.057
S      18.266  30.313  22.800
P      16.658  30.026  26.111
S      13.495  33.225  25.453
"""
if __name__ == "__main__":

    A = np.array( [21.436,  33.390,  24.057] , dtype=np.float128)
    B = np.array( [18.266,  30.313,  22.800] , dtype=np.float128)
    C = np.array( [16.658,  30.026,  26.111] , dtype=np.float128)
    D = np.array( [13.495,  33.225,  25.453] , dtype=np.float128)

    print ('Position of A = ', A)
    print ('Position of B = ', B)
    print ('Position of C = ', C)
    print ('Position of D = ', D)
    print ('')

    bond_AB = bond_length(A, B)
    bond_BC = bond_length(B, C)
    bond_CD = bond_length(C, D)
    print('bond(AB) = ', bond_AB)
    print('bond(BC) = ', bond_BC)
    print('bond(CD) = ', bond_CD)
    print ('')

    angl_ABC = bond_angle(A, B, C)
    angl_BCD = bond_angle(B, C, D)
    print('angle(ABC) = ', angl_ABC, ' = ', angl_ABC/np.pi*180.0)
    print('angle(BCD) = ', angl_BCD, ' = ', angl_BCD/np.pi*180.0)
    print ('')

    dihd_ABCD = torsion(A, B, C, D)
    print('dihedral(ABCD) = ', dihd_ABCD, ' = ', dihd_ABCD/np.pi*180.0)
    print ('')


    print('########################################')
    print('Calculate the coordinates of D by NeRF')

    D_test = NeRF(A, B, C, bond_CD, angl_BCD, dihd_ABCD)
    dist = np.linalg.norm(D - D_test)

    print('New Position of D = ', D_test)
    print('(Distance to the original D = ', dist, ')')
    print ('')

    new_bond_CD = bond_length(C, D_test)
    new_angl_BCD = bond_angle(B, C, D_test)
    new_dihd_ABCD = torsion(A, B, C, D_test)
    print('New bond(CD) = ', new_bond_CD)
    print('       Error = ', new_bond_CD - bond_CD)
    print('New angle(BCD) = ', new_angl_BCD, ' = ', new_angl_BCD/np.pi*180.0)
    print('         Error = ', new_angl_BCD/np.pi*180.0 - angl_BCD/np.pi*180.0 , ' degree')
    print('New dihedral(ABCD) = ', new_dihd_ABCD, ' = ', new_dihd_ABCD/np.pi*180.0)
    print('             Error = ', new_dihd_ABCD/np.pi*180.0 - dihd_ABCD/np.pi*180.0, ' degree')

