#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Convert VTF format to PSF + DCD files.

Usage:
    python vtf2psf_dcd.py -i input.vtf -o output_prefix

Produces output_prefix.psf and output_prefix.dcd
'''

import argparse
from lop.file_io.psf import PsfFile
from lop.file_io.dcd import DcdFile, DcdHeader
from lop.elements.psf import Psf, Atom


def parse_vtf_header(f):
    """Parse the VTF header (everything before the first 'timestep' line).
    Returns (atoms_info, bonds, pbc, header_end_pos).
    atoms_info: list of (name, type_str) indexed by atom index (0-based)
    bonds: list of (a, b) 0-based
    pbc: [x, y, z] or None
    """
    atoms_info = {}  # {index: (name, type_str)}
    bonds = []
    pbc = None

    while True:
        pos = f.tell()
        line = f.readline()
        if not line:
            break
        line = line.strip()

        if line.startswith('timestep'):
            f.seek(pos)
            break

        if not line:
            continue

        if line.startswith('pbc'):
            parts = line.split()
            pbc = [float(parts[1]), float(parts[2]), float(parts[3])]

        elif line.startswith('atom'):
            # e.g. "atom 0:4999 radius 0.5 name O type 0"
            parts = line.split()
            range_str = parts[1]
            # Parse keyword arguments
            kw = {}
            i = 2
            while i < len(parts) - 1:
                kw[parts[i]] = parts[i + 1]
                i += 2
            name = kw.get('name', 'X')
            type_str = kw.get('type', '0')
            # Expand range
            if ':' in range_str:
                start, end = range_str.split(':')
                for idx in range(int(start), int(end) + 1):
                    atoms_info[idx] = (name, type_str)
            else:
                atoms_info[int(range_str)] = (name, type_str)

        elif line.startswith('bond'):
            # e.g. "bond 1:0"
            parts = line.split()
            pair = parts[1].split(':')
            a, b = int(pair[0]), int(pair[1])
            bonds.append((a, b))

    n_atoms = max(atoms_info.keys()) + 1 if atoms_info else 0
    atoms_list = []
    for i in range(n_atoms):
        atoms_list.append(atoms_info.get(i, ('X', '0')))

    return atoms_list, bonds, pbc


def build_psf(atoms_list, bonds):
    """Build a Psf object from parsed VTF data."""
    psf = Psf()
    for i, (name, type_str) in enumerate(atoms_list):
        atom = Atom()
        atom.atom_id = i + 1  # 1-indexed
        atom.seg_name = 'A'
        atom.res_id = i + 1
        atom.res_name = 'POL'
        atom.atom_name = name
        atom.atom_type = name
        atom.charge = 0.0
        atom.mass = 1.0
        atom.unused = 0
        psf.atoms.append(atom)

    # Convert bonds to 1-indexed
    for a, b in bonds:
        psf.bonds.append((a + 1, b + 1))

    return psf


def main():
    parser = argparse.ArgumentParser(description='Convert VTF to PSF + DCD')
    parser.add_argument('-i', required=True, help='Input VTF file')
    parser.add_argument('-o', required=True, help='Output prefix (produces PREFIX.psf and PREFIX.dcd)')
    args = parser.parse_args()

    # Parse VTF
    with open(args.i, 'r') as f:
        atoms_list, bonds, pbc = parse_vtf_header(f)
        n_atoms = len(atoms_list)
        print(f'Number of atoms: {n_atoms}')
        print(f'Number of bonds: {len(bonds)}')
        if pbc:
            print(f'PBC: {pbc}')

        # Write PSF
        psf = build_psf(atoms_list, bonds)
        psf_file = PsfFile(args.o + '.psf')
        psf_file.open_to_write()
        psf_file.write_all(psf)
        psf_file.close()
        print(f'Wrote {args.o}.psf')

        # Count frames and read coordinates
        frames = []
        while True:
            line = f.readline()
            if not line:
                break
            line = line.strip()
            if line.startswith('timestep'):
                coord = []
                for _ in range(n_atoms):
                    cline = f.readline().split()
                    coord.append([float(cline[0]), float(cline[1]), float(cline[2])])
                frames.append(coord)

    n_frames = len(frames)
    print(f'Number of frames: {n_frames}')

    # Write DCD
    dcd = DcdFile(args.o + '.dcd')
    dcd.open_to_write()

    header = DcdHeader()
    header.nmp_real = n_atoms
    header.nset = n_frames
    header.istart = 0
    header.nstep_save = 1
    header.nstep = n_frames
    header.nunit_real = 0
    header.delta = 1.0
    header.tempk = 0.0
    header.lunit2mp = []
    header.title = [b'Created by vtf2psf_dcd.py'.ljust(80), b''.ljust(80)]
    header.format = 'charmm'
    if pbc:
        header.with_unit_cell = True
        header.unit_cell_xyz = pbc
        header.unit_cell_abc = [90.0, 90.0, 90.0]
    # block1: 21 elements (4s + 20i), index 20 = 24 for VMD compatibility
    header.block1 = [b'CORD'] + [0] * 19 + [24]

    dcd.set_header(header)
    dcd.write_header()

    for coord in frames:
        dcd.write_onestep(coord)

    dcd.close()
    print(f'Wrote {args.o}.dcd')


if __name__ == '__main__':
    main()
