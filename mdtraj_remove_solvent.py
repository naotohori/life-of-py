#!/usr/bin/env python

import sys

import mdtraj as md
import numpy as np

# Target solvent atomnames.
# Using resname does not work becuase mdtraj automatically converts WAT to HOH
solvent_atomnames = ['WAT', 'ETH']

if len(sys.argv) != 4:
    print('Usage: SCRIPT (input topology, e.g. psf) (input dcd) (output dcd)')
    sys.exit(2)

filepath_top = sys.argv[1]
filepath_dcd_in = sys.argv[2]
filepath_dcd_out = sys.argv[3]

traj = md.load(filepath_dcd_in, top=filepath_top)

# Get indeces
remove_atom_ids = [atom.index for atom in traj.topology.atoms if atom.name in solvent_atomnames]

# Indeces to leave
keep_atom_ids = np.setdiff1d(np.arange(traj.n_atoms), remove_atom_ids)

new_traj = traj.atom_slice(keep_atom_ids)

new_traj.save_dcd(filepath_dcd_out)

#new_traj.save_pdb('solute.pdb')

