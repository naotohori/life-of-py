#!/usr/bin/env python

import sys
import numpy as np
import argparse

from lop.file_io.xyz import XyzFile

def rotate_to_lay_at_origin(coords):
    # move the Center of Mass to the origin
    com = coords.mean(axis=0)
    coords = coords - com

    # gyration tensor
    G = (coords.T @ coords) / len(coords)

    eigvals, eigvecs = np.linalg.eigh(G)

    # Sort the eigenvectors in descending order
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]

    coords_rotated = coords @ eigvecs

    return coords_rotated


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Read XYZ file and lay on the stage at z = 0.',
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--stage-margin', dest='stage_margin', default=10.0,
                        action='store', type=float, 
                        help='The minimum distance between the stage (z=0) and the bottom of the molecule.')

    parser.add_argument('xyz_in', help='Input xyz file')
    parser.add_argument('xyz_out', help='Output xyz file')

    args = parser.parse_args()

    # Read xyz
    seq, coords = XyzFile(args.xyz_in, openmode='r').read(np_array=True, close=True)

    # Lay at origin
    coords_rot = rotate_to_lay_at_origin(coords)

    # Shift up the Z coordinates
    zmin = min(coords_rot[:,2])
    coords_rot[:,2] -= zmin - args.stage_margin

    # Write xyz
    XyzFile(args.xyz_out, openmode='w').write(seq, coords_rot, close=True)
