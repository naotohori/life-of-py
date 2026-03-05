#!/usr/bin/env python3

import sys

def compute_net_charge(mol2_file):
    total_charge = 0.0
    in_atom_section = False

    with open(mol2_file, "r") as f:
        for line in f:
            line = line.strip()

            if line.startswith("@<TRIPOS>ATOM"):
                in_atom_section = True
                continue

            if line.startswith("@<TRIPOS>") and in_atom_section:
                break

            if in_atom_section and line:
                parts = line.split()
                try:
                    charge = float(parts[-1])
                    total_charge += charge
                except (ValueError, IndexError):
                    pass

    return total_charge


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: SCRIPT molecule.mol2")
        sys.exit(1)

    mol2_file = sys.argv[1]
    net_charge = compute_net_charge(mol2_file)

    print(f"Net charge: {net_charge:.6f}")
