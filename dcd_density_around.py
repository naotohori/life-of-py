#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compute a 3D density map of molecule B around molecule A.

For each (frame, mol-A copy) pair:
  1. Unwrap mol-A beads (sequential MAXD algorithm, PBC).
  2. Superimpose mol-A onto a native reference via CalcROT → 4×4 matrix.
  3. Find mol-B chains that may have atoms inside the cubic box [-R, R]^3:
     conservative filter — skip a chain only if its centre is farther than
     R + chain_max_radius from the mol-A centre (minimum-image convention).
     Only atoms that actually fall inside [-R, R]^3 are binned.
  4. Apply the same 4×4 transform to nearby mol-B atoms.
  5. Bin transformed mol-B atom positions onto a cubic grid [-R, R]^3.

Output: OpenDX density file (atoms/Å³).
        Optional sequential PDB movie file (one MODEL per mol-A copy per frame).

Author: Naoto Hori, assisted by Claude Code
"""

import argparse
import sys
import math
import numpy as np
from numpy import zeros, float64

from lop.file_io.dcd import DcdFile
from lop.file_io.pdb import PdbFile

# ------------------------------------------------------------
# CalcROT import
# https://github.com/naotohori/fQCP
# ------------------------------------------------------------
from fQCP.CalcROT import calcrotation

# ------------------------------------------------------------
# Constants
# ------------------------------------------------------------
MAXD = 25.0  # Å — max jump between adjacent CG beads for unwrap


# ============================================================
# Helper functions
# ============================================================

def chain_atom_count(chain):
    return sum(len(r.atoms) for r in chain.residues)


def build_atom_slices(chains):
    """Return list of slice objects, one per chain, indexing into the flat atom array."""
    slices = []
    start = 0
    for chain in chains:
        n = chain_atom_count(chain)
        slices.append(slice(start, start + n))
        start += n
    return slices


def build_atom_meta(chains):
    """
    Return a list (one entry per chain) of lists of atom metadata tuples:
      (name, res_name, chain_id, res_seq, ins_code)
    Used for writing PDB ATOM lines with original identifiers.
    """
    meta = []
    for chain in chains:
        chain_meta = []
        for residue in chain.residues:
            for atom in residue.atoms:
                chain_meta.append((atom.name, atom.res_name, atom.chain_id,
                                   atom.res_seq, atom.ins_code))
        meta.append(chain_meta)
    return meta


def get_chain_id(chain):
    """Return the PDB chain_id of the first atom in a chain."""
    return chain.residues[0].atoms[0].chain_id if chain.residues else '?'


def chain_max_radius(chain):
    """
    Return the maximum distance from the geometric centre to any atom in the chain,
    computed from the static PDB coordinates.  Used as a conservative spatial filter:
    a chain can only have atoms inside a box of half-size R if its centre is within
    R + max_radius of the box centre.
    """
    xyz = np.array([atom.xyz.get_as_tuple()
                    for residue in chain.residues
                    for atom in residue.atoms], dtype=float64)
    centre = xyz.mean(axis=0)
    return float(np.sqrt(((xyz - centre) ** 2).sum(axis=1)).max())


def write_pdb_atom(fout, serial, name, res_name, chain_id, res_seq, ins_code, xyz,
                   occupancy=1.00, temp_factor=0.00):
    """Write one PDB ATOM line.

    ins_code is stored by the PDB parser as line[26:30] (4 chars: the insertion
    code byte plus the 3 mandatory blank columns 28-30), so {:4s} is correct here.
    """
    x, y, z = xyz
    fout.write(
        'ATOM  {:5d} {:4s}{:1s}{:3s} {:1s}{:4d}{:4s}'
        '{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}\n'.format(
            serial, name, ' ', res_name, chain_id, res_seq, ins_code,
            x, y, z, occupancy, temp_factor))


def write_movie_model(fout, model_num, frame_num_1based, mol_a_chain_id, mol_a_chain_idx,
                      mol_a_meta, mol_a_coords_tr,
                      mol_b_meta_list, mol_b_coords_tr_list):
    """
    Write one MODEL block to the sequential PDB movie file.

    Parameters
    ----------
    fout                 : open file handle
    model_num            : MODEL serial number (1-based, increments across all snapshots)
    frame_num_1based     : original DCD frame number (1-based)
    mol_a_chain_id       : PDB chain_id of this mol-A copy
    mol_a_chain_idx      : 1-based chain index of this mol-A copy in the system PDB
    mol_a_meta           : list of (name, res_name, chain_id, res_seq, ins_code) for mol-A atoms
    mol_a_coords_tr      : (N_a, 3) transformed mol-A coordinates
    mol_b_meta_list      : list of per-mol-B-chain meta lists (those inside the box)
    mol_b_coords_tr_list : list of (N_b, 3) transformed mol-B coordinate arrays
    """
    fout.write('MODEL     {:4d}\n'.format(model_num))
    fout.write('REMARK    DCD_FRAME {:d}  MOL_A_CHAIN_IDX {:d}  MOL_A_CHAIN_ID {:s}\n'.format(
        frame_num_1based, mol_a_chain_idx, mol_a_chain_id))

    serial = 1

    # mol-A atoms
    for i, (name, res_name, chain_id, res_seq, ins_code) in enumerate(mol_a_meta):
        write_pdb_atom(fout, serial, name, res_name, chain_id, res_seq, ins_code,
                       mol_a_coords_tr[i])
        serial += 1
    fout.write('TER\n')

    # mol-B atoms
    for meta, coords in zip(mol_b_meta_list, mol_b_coords_tr_list):
        for i, (name, res_name, chain_id, res_seq, ins_code) in enumerate(meta):
            write_pdb_atom(fout, serial, name, res_name, chain_id, res_seq, ins_code,
                           coords[i])
            serial += 1
        fout.write('TER\n')

    fout.write('ENDMDL\n')


def unwrap_molecule(coords, box):
    """
    Unwrap a single molecule using the sequential MAXD algorithm.
    coords : (N, 3) float array
    box    : [Lx, Ly, Lz]
    Returns unwrapped (N, 3) float64 array.
    """
    box = np.asarray(box, dtype=float64)
    n = len(coords)
    if n == 0:
        return coords
    result = coords.astype(float64)
    add = np.zeros(3, dtype=float64)
    pre = result[0].copy()
    for idx in range(1, n):
        xyz  = result[idx] + add          # new array; does not alias result[idx]
        diff = xyz - pre
        corr = np.where(diff > MAXD, -box, np.where(diff < -MAXD, box, 0.0))
        add += corr
        xyz += corr
        pre  = xyz                        # safe: next iteration creates a fresh xyz
        result[idx] = xyz
    return result


def nearest_image_shift(center_b, center_a, box):
    """
    Return the shift vector (Å) that brings center_b closest to center_a
    under the minimum-image convention.
    """
    box   = np.asarray(box, dtype=float64)
    delta = center_b - center_a
    return -np.round(delta / box) * box


def apply_transform(mat, coords):
    """
    Apply 4×4 homogeneous transform mat to coords.
    coords : (N, 3)
    Returns (N, 3) float64.
    Equivalent to the augmented-matrix form but avoids allocating a (4, N) array.
    """
    return coords @ mat[:3, :3].T + mat[:3, 3]


def write_opendx(filename, grid, origin, grid_size):
    """
    Write a 3D numpy array as an OpenDX file.
    grid      : (NX, NY, NZ) array
    origin    : (OX, OY, OZ) — corner of voxel [0,0,0]
    grid_size : voxel spacing in Å
    """
    nx, ny, nz = grid.shape
    total = nx * ny * nz
    flat = grid.flatten(order='C')

    with open(filename, 'w') as fout:
        fout.write('object 1 class gridpositions counts {:d} {:d} {:d}\n'.format(nx, ny, nz))
        fout.write('origin {:.6f} {:.6f} {:.6f}\n'.format(*origin))
        fout.write('delta {:.6f} 0 0\n'.format(grid_size))
        fout.write('delta 0 {:.6f} 0\n'.format(grid_size))
        fout.write('delta 0 0 {:.6f}\n'.format(grid_size))
        fout.write('object 2 class gridconnections counts {:d} {:d} {:d}\n'.format(nx, ny, nz))
        fout.write('object 3 class array type double rank 0 items {:d} data follows\n'.format(total))
        # Write 3 values per line using numpy for bulk formatting (faster than a Python loop).
        # Pad to a multiple of 3 with zeros; DX readers use the item count from the header
        # and ignore any trailing values, so padding is safe.
        remainder = total % 3
        if remainder:
            padded = np.concatenate([flat, np.zeros(3 - remainder, dtype=float64)])
        else:
            padded = flat
        np.savetxt(fout, padded.reshape(-1, 3), fmt='%.6e', delimiter='  ')
        fout.write('attribute "dep" string "positions"\n')
        fout.write('object "density" class field\n')
        fout.write('component "positions" value 1\n')
        fout.write('component "connections" value 2\n')
        fout.write('component "data" value 3\n')


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description=(
            'Compute a 3D density map of molecule B around molecule A '
            'from a PBC-wrapped DCD trajectory. '
            'Mol-A copies are superimposed onto a native reference; '
            'nearby mol-B chains are accumulated on a cubic grid.'))
    parser.add_argument('pdb',    help='System PDB file')
    parser.add_argument('dcd',    help='DCD trajectory file (PBC-wrapped)')
    parser.add_argument('native', help='Native mol-A PDB (single copy, superposition reference)')
    parser.add_argument('--mol_a', type=int, nargs=2, required=True,
                        metavar=('FIRST', 'LAST'),
                        help='1-based inclusive chain range for molecule A '
                             '(reference for superposition), e.g. --mol_a 1 40')
    parser.add_argument('--mol_b', type=int, nargs=2, default=None,
                        metavar=('FIRST', 'LAST'),
                        help='1-based inclusive chain range for molecule B '
                             '(density subject). Default: all chains not in mol-A.')
    parser.add_argument('--range',     type=float, required=True, dest='half_box',
                        help='Half-box size in Å; cubic analysis region [-R,R]^3 '
                             'centred on the superimposed mol-A')
    parser.add_argument('--grid_size', type=float, required=True,
                        help='Voxel spacing in Å')
    parser.add_argument('--output',    required=True,
                        help='Output OpenDX density filename')
    parser.add_argument('--frames',    type=int, nargs=2, metavar=('START', 'END'),
                        default=None,
                        help='1-based inclusive frame range (default: all frames)')
    parser.add_argument('--movie',     default=None, metavar='FILE',
                        help='Output sequential PDB movie file '
                             '(one MODEL per mol-A copy per frame)')
    args = parser.parse_args()

    R     = args.half_box
    dx    = args.grid_size
    ngrid = int(math.floor(2.0 * R / dx))
    origin = np.array([-R, -R, -R], dtype=float64)

    print('Grid: {:d}^3 voxels, spacing {:.2f} Å, half-box {:.2f} Å'.format(ngrid, dx, R))

    # ----------------------------------------------------------
    # Read system PDB → per-chain slices and metadata
    # ----------------------------------------------------------
    sys_pdb = PdbFile(args.pdb)
    sys_pdb.open_to_read()
    sys_chains = sys_pdb.read_all()
    sys_pdb.close()

    n_chains = len(sys_chains)

    # Convert 1-based inclusive ranges to 0-based index lists
    mol_a_first, mol_a_last = args.mol_a
    if not (1 <= mol_a_first <= mol_a_last <= n_chains):
        print('ERROR: --mol_a {:d} {:d} is out of range (system has {:d} chains)'.format(
            mol_a_first, mol_a_last, n_chains))
        sys.exit(1)
    mol_a_indices = list(range(mol_a_first - 1, mol_a_last))  # 0-based

    if args.mol_b is not None:
        mol_b_first, mol_b_last = args.mol_b
        if not (1 <= mol_b_first <= mol_b_last <= n_chains):
            print('ERROR: --mol_b {:d} {:d} is out of range (system has {:d} chains)'.format(
                mol_b_first, mol_b_last, n_chains))
            sys.exit(1)
        mol_b_indices = list(range(mol_b_first - 1, mol_b_last))  # 0-based
    else:
        mol_a_set = set(mol_a_indices)
        mol_b_indices = [i for i in range(n_chains) if i not in mol_a_set]

    num_mol_a = len(mol_a_indices)
    num_mol_b = len(mol_b_indices)

    print('System PDB: {:d} chains total  |  mol-A: {:d} chains (indices {:d}–{:d})  '
          '|  mol-B: {:d} chains'.format(
              n_chains, num_mol_a, mol_a_indices[0] + 1, mol_a_indices[-1] + 1, num_mol_b))

    all_slices = build_atom_slices(sys_chains)
    all_meta   = build_atom_meta(sys_chains)

    mol_a_slices   = [all_slices[i] for i in mol_a_indices]
    mol_a_meta_all = [all_meta[i]   for i in mol_a_indices]
    mol_a_chain_ids = [get_chain_id(sys_chains[i]) for i in mol_a_indices]

    mol_b_slices   = [all_slices[i] for i in mol_b_indices]
    mol_b_meta_all = [all_meta[i]   for i in mol_b_indices]

    # Conservative spatial filter: a mol-B chain can contribute atoms inside the
    # analysis box only if its centre is within R + mol_b_global_radius of the mol-A
    # centre.  Because mol-B is flexible, we use the single largest chain_max_radius
    # across all mol-B chains (from the static PDB) as a universal upper bound, rather
    # than per-chain radii that may underestimate the true extent during simulation.
    mol_b_max_radii     = np.array([chain_max_radius(sys_chains[i]) for i in mol_b_indices])
    mol_b_global_radius = float(mol_b_max_radii.max())
    print('Mol-B chain radii (static PDB): min={:.1f} Å  max={:.1f} Å  '
          '(using {:.1f} Å as universal filter threshold)'.format(
              mol_b_max_radii.min(), mol_b_max_radii.max(), mol_b_global_radius))

    # ----------------------------------------------------------
    # Read native mol-A → reference coords (3, N_a) Fortran order
    # ----------------------------------------------------------
    nat_pdb = PdbFile(args.native)
    nat_pdb.open_to_read()
    nat_chains = nat_pdb.read_all()
    nat_pdb.close()

    ref_list = []
    for chain in nat_chains:
        for residue in chain.residues:
            for atom in residue.atoms:
                ref_list.append(atom.xyz.get_as_tuple())
    ref_np   = np.array(ref_list, dtype=float64)   # (N_a, 3)
    ref_mol_a_F = np.asfortranarray(ref_np.T)       # (3, N_a) Fortran order

    center_native = ref_np.mean(axis=0)
    N_mol_a_atoms = ref_np.shape[0]
    print('Native mol-A: {:d} atoms, center ({:.3f}, {:.3f}, {:.3f})'.format(
        N_mol_a_atoms, *center_native))

    # Verify all mol-A copies have the same atom count as the native
    for rank, (i, sl) in enumerate(zip(mol_a_indices, mol_a_slices)):
        n = sl.stop - sl.start
        if n != N_mol_a_atoms:
            print('ERROR: mol-A chain {:d} (system chain {:d}) has {:d} atoms, '
                  'native has {:d}'.format(rank + 1, i + 1, n, N_mol_a_atoms))
            sys.exit(1)

    # ----------------------------------------------------------
    # Open DCD
    # ----------------------------------------------------------
    dcd = DcdFile(args.dcd)
    dcd.open_to_read()
    dcd.read_header()
    header = dcd.get_header()
    nframes_total = dcd.count_frame()
    print('DCD frames: {:d}, atoms: {:d}'.format(nframes_total, header.nmp_real))

    # Frame range (1-based → 0-based internally)
    if args.frames is not None:
        frame_start = args.frames[0] - 1
        frame_end   = args.frames[1] - 1
    else:
        frame_start = 0
        frame_end   = nframes_total - 1

    # ----------------------------------------------------------
    # Density grid
    # ----------------------------------------------------------
    grid          = zeros((ngrid, ngrid, ngrid), dtype=float64)
    total_samples = 0

    # ----------------------------------------------------------
    # Open movie file (optional)
    # ----------------------------------------------------------
    movie_fout = open(args.movie, 'w') if args.movie else None
    model_num  = 0

    # ----------------------------------------------------------
    # Per-frame loop
    # ----------------------------------------------------------
    dcd.rewind()
    if frame_start > 0:
        dcd.skip(frame_start)

    for iframe in range(frame_start, frame_end + 1):
        if not dcd.has_more_data():
            print('Warning: DCD ended at frame {:d}'.format(iframe + 1))
            break

        coords = dcd.read_onestep_np()      # (N_total, 3)
        box    = dcd._header.unit_cell_xyz  # [Lx, Ly, Lz]
        box_np = np.asarray(box, dtype=float64)

        frame_rmsd_list = []

        for rank_a, sl_a in enumerate(mol_a_slices):

            # ---- Extract and unwrap mol-A copy ----
            mol_a_raw   = coords[sl_a].astype(float64)
            mol_a_unwrp = unwrap_molecule(mol_a_raw, box_np)
            mol_a_center = mol_a_unwrp.mean(axis=0)

            # ---- Superimpose onto native ----
            mol_a_F = np.asfortranarray(mol_a_unwrp.T)          # (3, N_a)
            rmsd, mat = calcrotation(ref_mol_a_F, mol_a_F)      # mat: 4×4
            frame_rmsd_list.append(rmsd)

            # Transform mol-A itself (needed for movie)
            mol_a_tr = apply_transform(mat, mol_a_unwrp) if movie_fout else None

            # ---- Find and process nearby mol-B chains ----
            mol_b_meta_in_box   = []
            mol_b_coords_in_box = []
            mol_b_tr_batch      = []   # accumulate transformed coords for bulk binning

            for rank_b, sl_b in enumerate(mol_b_slices):
                mol_b_raw        = coords[sl_b].astype(float64)
                mol_b_center_raw = mol_b_raw.mean(axis=0)

                # Minimum-image shift to bring mol-B center near mol-A center
                shift        = nearest_image_shift(mol_b_center_raw, mol_a_center, box_np)
                mol_b_center = mol_b_center_raw + shift

                # Conservative filter: skip only if the chain is entirely outside the
                # box, i.e. even its nearest atom cannot reach [-R, R]^3.
                # A chain with centre distance d and radius r has no atoms in the box
                # when d - r > R in any dimension.
                if np.any(np.abs(mol_b_center - mol_a_center) > R + mol_b_global_radius):
                    continue

                # Apply same image shift to all mol-B atoms, then unwrap
                mol_b_coords = unwrap_molecule(mol_b_raw + shift, box_np)

                # Apply the 4×4 superposition transform
                mol_b_tr = apply_transform(mat, mol_b_coords)   # (N_b, 3)

                mol_b_tr_batch.append(mol_b_tr)

                if movie_fout:
                    mol_b_meta_in_box.append(mol_b_meta_all[rank_b])
                    mol_b_coords_in_box.append(mol_b_tr)

            # Vectorised binning — all nearby mol-B atoms in a single numpy pass
            if mol_b_tr_batch:
                all_tr = np.vstack(mol_b_tr_batch)
                ixyz   = np.floor((all_tr - origin) / dx).astype(np.intp)
                mask   = np.all((ixyz >= 0) & (ixyz < ngrid), axis=1)
                valid  = ixyz[mask]
                if len(valid):
                    flat_idx = (valid[:, 0] * ngrid + valid[:, 1]) * ngrid + valid[:, 2]
                    grid += np.bincount(flat_idx,
                                        minlength=ngrid * ngrid * ngrid
                                        ).reshape(ngrid, ngrid, ngrid)

            # ---- Write movie snapshot ----
            if movie_fout:
                model_num += 1
                write_movie_model(
                    movie_fout, model_num, iframe + 1,
                    mol_a_chain_ids[rank_a], mol_a_indices[rank_a] + 1,
                    mol_a_meta_all[rank_a], mol_a_tr,
                    mol_b_meta_in_box, mol_b_coords_in_box)

            total_samples += 1

        avg_rmsd = np.mean(frame_rmsd_list) if frame_rmsd_list else float('nan')
        print('Frame {:5d}  avg_RMSD={:.3f} Å'.format(iframe + 1, avg_rmsd))

    dcd.close()

    if movie_fout:
        movie_fout.close()
        print('Movie written to {:s}  ({:d} models)'.format(args.movie, model_num))

    print('\nTotal (frame × mol-A copy) samples: {:d}'.format(total_samples))

    # ----------------------------------------------------------
    # Normalize: atoms/Å³
    # ----------------------------------------------------------
    if total_samples > 0:
        density = grid / (total_samples * dx**3)
    else:
        density = grid

    # ----------------------------------------------------------
    # Write OpenDX
    # ----------------------------------------------------------
    write_opendx(args.output, density, tuple(origin), dx)
    print('Density written to', args.output)
    integral = density.sum() * dx**3
    print('Integral over grid (avg mol-B atoms in box per sample): {:.3f}'.format(integral))


if __name__ == '__main__':
    main()
