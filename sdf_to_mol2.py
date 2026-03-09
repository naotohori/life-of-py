#!/usr/bin/env python3

import sys
import os
from openbabel import pybel


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} output.mol2 input1.sdf [input2.sdf ...]", file=sys.stderr)
        sys.exit(1)

    output_mol2 = sys.argv[1]
    input_sdfs = sys.argv[2:]

    # Check the input file
    for sdf in input_sdfs:
        if not os.path.isfile(sdf):
            print(f"Error: input file not found: {sdf}", file=sys.stderr)
            sys.exit(1)

    count = 0

    # Load sdf files and write to mol2, keeping the order
    with open(output_mol2, "w") as out_f:
        for sdf in input_sdfs:
            try:
                for mol in pybel.readfile("sdf", sdf):
                    out_f.write(mol.write("mol2"))
                    count += 1
            except Exception as e:
                print(f"Error while processing {sdf}: {e}", file=sys.stderr)
                sys.exit(1)

    print(f"Wrote {count} structures to {output_mol2}")


if __name__ == "__main__":
    main()
