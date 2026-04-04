#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Generate MOL (V2000) file from PSF and XYZ files.

Usage:
    python psf_xyz2mol.py -p input.psf -x input.xyz -o output.mol
'''

import sys
import argparse

from lop.file_io.psf import PsfFile
from lop.file_io.xyz import XyzFile


def main():
    parser = argparse.ArgumentParser(description='Generate MOL (V2000) from PSF and XYZ')
    parser.add_argument('-p', '--psf', required=True, help='Input PSF file')
    parser.add_argument('-x', '--xyz', required=True, help='Input XYZ file')
    parser.add_argument('-o', '--out', required=True, help='Output MOL file')
    args = parser.parse_args()

    # Read PSF
    psf_file = PsfFile(args.psf)
    psf_file.open_to_read()
    psf = psf_file.read_all()
    psf_file.close()

    # Read XYZ
    seq, coords = XyzFile(args.xyz, openmode='r').read(close=True)

    n_atoms = len(psf.atoms)
    n_bonds = len(psf.bonds)

    if n_atoms != len(coords):
        print(f'Error: PSF has {n_atoms} atoms but XYZ has {len(coords)} atoms')
        sys.exit(2)

    atom_names = [atom.atom_name.strip() for atom in psf.atoms]

    # Write MOL V2000
    with open(args.out, 'w') as f:
        # Header block (3 lines)
        f.write('\n')  # molecule name
        f.write('\n')  # program/timestamp
        f.write('\n')  # comment

        # Counts line
        f.write('%3d%3d  0  0  0  0  0  0  0  0999 V2000\n' % (n_atoms, n_bonds))

        # Atom block
        for i in range(n_atoms):
            x, y, z = coords[i]
            f.write('%10.4f%10.4f%10.4f %-3s 0  0  0  0  0  0  0  0  0  0  0  0\n'
                    % (x, y, z, atom_names[i]))

        # Bond block
        for a, b in psf.bonds:
            f.write('%3d%3d%3d%3d%3d%3d%3d\n' % (a, b, 1, 0, 0, 0, 0))

        # End
        f.write('M  END\n')

    print(f'Wrote {args.out} ({n_atoms} atoms, {n_bonds} bonds)')


if __name__ == '__main__':
    main()
