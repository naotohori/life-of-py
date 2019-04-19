#!/usr/bin/env python

from cafysis.file_io.dcd import DcdFile
from cafysis.file_io.pdb import PdbFile
import py_bestfit
from numpy import zeros, asarray, float64
import sys

if len(sys.argv) < 5:
    print(('Usage: % SCRIPT [PDB filename] [DCD filename]'
            +'[serial ID begin] [serial ID end]'))
    sys.exit(2)
    
filename_pdb = sys.argv[1]
filename_dcd = sys.argv[2]

id_begin = []
id_end = []
for iarg in range(3, len(sys.argv)-1, 2) :
    id_begin.append(int(sys.argv[iarg]))
    id_end.append(int(sys.argv[iarg+1]))

# Coord1
pdb = PdbFile(filename_pdb)
pdb.open_to_read()
ref_chains = pdb.read_all()
pdb.close()

num_atom = 0
for chain in ref_chains :
    for residue in chain.residues :
        num_atom += len(residue.atoms)
    
ref = zeros((3, num_atom), dtype=float64, order='F')

i = 0
for chain in ref_chains :
    for residue in chain.residues:
        for atom in residue.atoms :
            (ref[0][i], ref[1][i], ref[2][i]) = atom.xyz.get_as_tuple()
            i += 1

ref_idx = []
pre_idx = []

# all to all
for i in range(len(id_begin)) :
    for j in range(id_begin[i], id_end[i]+1) :
        ref_idx.append(j)
        pre_idx.append(j)
    
dcd = DcdFile(filename_dcd)
dcd.open_to_read()
dcd.read_header()

#dcd.show_header()
k = 0
while dcd.has_more_data() :
    k += 1
    pre_dcd = dcd.read_onestep()
    pre = asarray(zeros((3, num_atom), dtype=float64, order='F'))
    for i in range(3):
        for j in range(num_atom) :
            pre[i][j] = pre_dcd[j][i]
    
    (post, rmsd, ier, rot, center_ref, center_pre) = py_bestfit.bestfit(ref, pre,
                                                                        ref_idx, pre_idx)
    
    if ier == 0 :
        print('%5.2f' % rmsd)
    else :
        print('error in %i-th structure' % k)
        
dcd.close()

## F2PY
#bestfit - Function signature:
#  coord3,rmsd,ierr,r,xc1,xc2 = bestfit(coord1,coord2,list1,list2,[nat1,nat2,nat])
#Required arguments:
#  coord1 : input rank-2 array('d') with bounds (3,nat1)
#  coord2 : input rank-2 array('d') with bounds (3,nat2)
#  list1 : input rank-1 array('i') with bounds (nat)
#  list2 : input rank-1 array('i') with bounds (nat)
#Optional arguments:
#  nat1 := shape(coord1,1) input int
#  nat2 := shape(coord2,1) input int
#  nat := len(list1) input int
#Return objects:
#  coord3 : rank-2 array('d') with bounds (3,nat2)
#  rmsd : float
#  ierr : int
#  r : rank-2 array('d') with bounds (3,3)
#  xc1 : rank-1 array('d') with bounds (3)
#  xc2 : rank-1 array('d') with bounds (3)


#  4 !      bestfit.f            Version 1 5/19/1993      Patrice Koehl
#  5 !
#  6 !      This subroutine performs the best fit of two structures
#  7 !      (i.e. moves structure 2 such as to provide the best superposition
#  8 !      of the two structure)
#  9 !      The new structure 2 is called structure 3
# 10 !      It is based on the algorithm of Mclachlan
# 11 !
# 12 !      Input:
# 13 !            -coord1 (size 3xnat1) : vector containing the coordinates
# 14 !                              of all atoms of molecule 1
# 15 !                              arranged as x1,y1,z1,x2,y2,z2...
# 16 !                              Molecule 1 is considered the
# 17 !                              "reference" structure
# 18 !            -nat1                  :      number of atoms in molecule 1
# 19 !            -coord2 (size 3xnat2) : coordinates of all atoms of
# 20 !                              molecule 2
# 21 !                              Molecule 2 is considered the "test"
# 22 !                              structure; this is the one
# 23 !                              that will be moved 
# 24 !            -nat2                  : number of atoms in molecule 2
# 25 !            -nat                  : number of atoms picked in both
# 26 !                              structures for superposition
# 27 !            list1 (size nat)      :      position of the nat atoms for
# 28 !                              superposition in molecule 1
# 29 !            list2 (size nat)      :      position of the nat atoms for
# 30 !                              superposition in molecule 2
# 31 !      Output:
# 32 !            -coord3 (size 3xnat2) : coordinates of all atoms of molecule
# 33 !                              2 after superposition. coord2
# 34 !                              remains unchanged
# 35 !            -rmsd                  : coordinate RMS between the 2
# 36 !                              structure (caculated over the
# 37 !                              nat atoms picked)
# 38 !            -ierr                  : flag: 0 if computation of rmsd
# 39 !                              was OK, 1 otherwise
# 40 !
# 47 subroutine bestfit( coord1,  & ! 1[i ] ST1 structure
# 48                     nat1,    & ! 2[i ] # of ST1 atom
# 49                     coord2,  & ! 3[i ] ST2 structure
# 50                     nat2,    & ! 4[i ] # of ST2 atom
# 51                     nat,     & ! 5[i ] # of atom for calculation
# 52                     coord3,  & ! 6[ o] fitted structure
# 53                     list1,   & ! 7[i ] ST1 atom index for calculation
# 54                     list2,   & ! 8[i ] ST2 atom index for calculation
# 55                     rmsd,    & ! 9[ o] RMSD
# 56                     ierr,    & !10[ o] Error code
# 57                     r,       & !11[ o] Rotational matrix
# 58                     xc1,     & !12[ o] Center of ST1
# 59                     xc2 )      !13[ o] Center of ST2
# 60 
# 61    ! INPUT
# 62    integer, intent(in)  :: nat
# 63    integer, intent(in)  :: nat1
# 64    integer, intent(in)  :: nat2
# 65    integer, intent(in)  :: list1(nat)
# 66    integer, intent(in)  :: list2(nat)
# 67    real*8,  intent(in)  :: coord1(3,nat1)
# 68    real*8,  intent(in)  :: coord2(3,nat2)
# 69 
# 70    ! OUTPUT
# 71    integer, intent(out) :: ierr
# 72    real*8,  intent(out) :: coord3(3,nat2)
# 73    real*8,  intent(out) :: rmsd
# 74    real*8,  intent(out) :: r(3,3)  ! rotational matrix
# 75    real*8,  intent(out) :: xc1(3)
# 76    real*8,  intent(out) :: xc2(3)
