#!/usr/bin/env python3
"""
- Use RDKit to generate multiple possible conformers.
- Read a single-structure MOL2 (input reference coordinates)
- Add Hs (with coordinates)
- Iteratively generate conformers (ETKDGv3) WITHOUT RDKit pruning
- Optimise each conformer (MMFF94s if available, else UFF)
- Maintain a growing pool, sort by energy, and select RMSD-unique conformers greedily
- Ensure the FINAL number of output conformers equals -n (target), by generating more if needed
- Output: one SDF per conformer: <prefix>_<index>.sdf (index 0 = min energy)
- Output text: index, energy, RMSD_to_min, RMSD_to_input
  -o/--out-txt to set; default: <prefix>_energies.txt

RMSD:
- RMSD_to_min: aligned RMSD to the minimum-energy conformer among the FINAL selected set
- RMSD_to_input: aligned RMSD to the original input MOL2 coordinates (after adding H coords)

Notes:
- -n is the FINAL output count (RMSD-unique). Internally more conformers may be generated.
"""

import argparse
import os
from typing import List, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign


def mmff_available(mol: Chem.Mol) -> bool:
    try:
        return AllChem.MMFFHasAllMoleculeParams(mol)
    except Exception:
        return False


def optimize_and_energy(
    mol: Chem.Mol,
    cid: int,
    mmff_variant: str = "MMFF94s",
    max_iters: int = 200,
    use_mmff: bool = True,
) -> float:
    """Optimize one conformer and return FF energy (MMFF if possible, else UFF)."""
    e = float("inf")

    if use_mmff:
        try:
            AllChem.MMFFOptimizeMolecule(
                mol,
                mmffVariant=mmff_variant,
                maxIters=max_iters,
                confId=cid,
                ignoreInterfragInteractions=True,
            )
            ff = AllChem.MMFFGetMoleculeForceField(mol, confId=cid, mmffVariant=mmff_variant)
            if ff is not None:
                return float(ff.CalcEnergy())
        except Exception:
            pass

    # UFF fallback
    try:
        AllChem.UFFOptimizeMolecule(mol, confId=cid, maxIters=max_iters)
        ff = AllChem.UFFGetMoleculeForceField(mol, confId=cid)
        if ff is not None:
            e = float(ff.CalcEnergy())
    except Exception:
        e = float("inf")

    return e


def ensure_parent_dir(path_prefix: str) -> None:
    d = os.path.dirname(path_prefix)
    if d:
        os.makedirs(d, exist_ok=True)


def make_single_conformer_copy(mol: Chem.Mol, conf_id: int) -> Chem.Mol:
    """
    Make a molecule copy that contains exactly one conformer, with coordinates copied
    atom-by-atom. Robust against RDKit builds where conformer ownership causes errors.
    """
    m = Chem.Mol(mol)
    m.RemoveAllConformers()

    src_conf = mol.GetConformer(conf_id)
    new_conf = Chem.Conformer(m.GetNumAtoms())
    new_conf.SetId(0)
    for i in range(m.GetNumAtoms()):
        new_conf.SetAtomPosition(i, src_conf.GetAtomPosition(i))

    m.AddConformer(new_conf, assignId=True)
    return m


def write_one_conformer_sdf(
    mol: Chem.Mol,
    conf_id: int,
    out_path: str,
    name: str,
    energy: float,
) -> None:
    """Write one SDF containing one molecule with exactly one conformer."""
    m = make_single_conformer_copy(mol, conf_id)
    m.SetProp("_Name", name)
    m.SetProp("ENERGY", f"{energy:.6f}")

    w = Chem.SDWriter(out_path)
    try:
        w.write(m)
    finally:
        w.close()


def rmsd_aligned(prb_mol: Chem.Mol, prb_cid: int, ref_mol: Chem.Mol, ref_cid: int = 0) -> float:
    """
    Compute aligned RMSD (Å) between prb conformer and ref conformer without modifying originals.
    Uses AlignMol on a prb copy containing a single conformer.
    """
    prb = make_single_conformer_copy(prb_mol, prb_cid)
    try:
        return float(rdMolAlign.AlignMol(prb, ref_mol, prbCid=0, refCid=ref_cid))
    except TypeError:
        # older signature fallback
        return float(rdMolAlign.AlignMol(prb, ref_mol, 0, ref_cid))


def select_rmsd_unique_by_energy(
    mol: Chem.Mol,
    pool_sorted: List[Tuple[int, float]],
    target_n: int,
    rmsd_thresh: float,
) -> List[Tuple[int, float]]:
    """
    Greedy selection in ascending energy order.
    Keep a conformer if its RMSD to ALL already-kept conformers is >= rmsd_thresh.
    Returns up to target_n items.
    """
    selected: List[Tuple[int, float]] = []
    for cid, e in pool_sorted:
        if not selected:
            selected.append((cid, e))
            if len(selected) >= target_n:
                break
            continue

        keep = True
        for scid, _ in selected:
            r = rmsd_aligned(mol, cid, mol, ref_cid=scid)
            if r < rmsd_thresh:
                keep = False
                break
        if keep:
            selected.append((cid, e))
            if len(selected) >= target_n:
                break
    return selected


def main():
    ap = argparse.ArgumentParser(
        description="Generate RMSD-unique conformers (final count = -n) from input MOL2; write one SDF per conformer + energy/RMSD text."
    )
    ap.add_argument("-i", "--input", required=True, help="Input MOL2 (single structure).")
    ap.add_argument(
        "-p", "--prefix", required=True,
        help="Output prefix (may include folders). Output: <prefix>_<index>.sdf"
    )
    ap.add_argument(
        "-n", "--num-confs", type=int, default=100,
        help="FINAL number of RMSD-unique conformers to output."
    )
    ap.add_argument("--threads", type=int, default=0, help="Threads for embedding (ETKDG).")
    ap.add_argument(
        "--prune-rms", type=float, default=0.5,
        help="RMSD threshold (Å) used for UNIQUE selection (final set)."
    )
    ap.add_argument("--seed", type=int, default=1, help="Random seed.")
    ap.add_argument("--pad", type=int, default=4, help="Zero-pad width for index.")
    ap.add_argument("--mmff-variant", default="MMFF94s", help="MMFF variant: MMFF94 or MMFF94s.")
    ap.add_argument("--max-iters", type=int, default=200, help="Max optimization iterations.")
    ap.add_argument(
        "-o", "--out-txt", default=None,
        help="Output text filename. Default: <prefix>_energies.txt"
    )

    # New controls for oversampling / retries (defaults chosen to be reasonable)
    ap.add_argument(
        "--oversample", type=float, default=5.0,
        help="How many more conformers to generate per remaining needed (default: 5.0)."
    )
    ap.add_argument(
        "--min-batch", type=int, default=50,
        help="Minimum number of conformers to attempt per generation round (default: 50)."
    )
    ap.add_argument(
        "--max-rounds", type=int, default=20,
        help="Maximum generation rounds to try to reach -n unique conformers (default: 20)."
    )

    args = ap.parse_args()

    target_n = int(args.num_confs)
    if target_n < 1:
        raise SystemExit("ERROR: -n/--num-confs must be >= 1")

    mol = Chem.MolFromMol2File(args.input, sanitize=True, removeHs=False)
    if mol is None:
        raise SystemExit(f"ERROR: failed to read MOL2: {args.input}")

    # Add Hs and preserve/add coordinates (this defines input reference geometry)
    #mol = Chem.AddHs(mol, addCoords=True)
    if mol.GetNumConformers() < 1:
        raise SystemExit("ERROR: input molecule has no coordinates/conformer; cannot compute RMSD_to_input.")

    input_ref = make_single_conformer_copy(mol, conf_id=0)  # reference conformer id=0 in input_ref

    use_mmff = mmff_available(mol)

    # Pool: list of (confId, energy)
    pool: List[Tuple[int, float]] = []

    # Iteratively generate until we can select target_n unique conformers
    selected: List[Tuple[int, float]] = []

    for rnd in range(int(args.max_rounds)):
        need = target_n - len(selected)
        if need <= 0:
            break

        # Batch size: oversample * remaining need, but at least min-batch
        batch = max(int(round(need * float(args.oversample))), int(args.min_batch))

        params = AllChem.ETKDGv3()
        params.randomSeed = int(args.seed) + rnd  # vary per round deterministically
        params.numThreads = int(args.threads)
        # IMPORTANT: disable RDKit internal pruning so we control uniqueness ourselves
        params.pruneRmsThresh = -1.0

        new_cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=batch, params=params))
        if not new_cids:
            # If embedding fails completely, no point continuing
            break

        # Optimize and score only the new conformers
        for cid in new_cids:
            e = optimize_and_energy(
                mol, cid,
                mmff_variant=str(args.mmff_variant),
                max_iters=int(args.max_iters),
                use_mmff=use_mmff,
            )
            pool.append((cid, e))

        # Sort whole pool by energy and re-select from scratch (ensures global best among generated)
        pool_sorted = sorted(pool, key=lambda x: x[1])
        selected = select_rmsd_unique_by_energy(
            mol, pool_sorted, target_n=target_n, rmsd_thresh=float(args.prune_rms)
        )

        if len(selected) >= target_n:
            break

    if len(selected) < target_n:
        raise SystemExit(
            f"ERROR: Could only select {len(selected)} RMSD-unique conformers (target={target_n}).\n"
            f"Try lowering --prune-rms, increasing --max-rounds, or increasing --oversample."
        )

    # Now selected is energy-ascending and length == target_n
    min_cid = selected[0][0]

    ensure_parent_dir(args.prefix)
    txt_path = args.out_txt if args.out_txt is not None else f"{args.prefix}_energies.txt"
    txt_dir = os.path.dirname(txt_path)
    if txt_dir:
        os.makedirs(txt_dir, exist_ok=True)

    with open(txt_path, "w") as ftxt:
        ftxt.write("# index\tenergy\tRMSD_to_min\tRMSD_to_input\n")
        for idx, (cid, e) in enumerate(selected):
            out_sdf = f"{args.prefix}_{idx:0{int(args.pad)}d}.sdf"
            name = f"CONF_{idx:0{int(args.pad)}d}"
            write_one_conformer_sdf(mol, cid, out_sdf, name=name, energy=e)

            rmsd_min = 0.0 if cid == min_cid else rmsd_aligned(mol, cid, mol, ref_cid=min_cid)
            rmsd_inp = rmsd_aligned(mol, cid, input_ref, ref_cid=0)

            ftxt.write(f"{idx}\t{e:.6f}\t{rmsd_min:.4f}\t{rmsd_inp:.4f}\n")

    print(f"Wrote {len(selected)} conformers to: {args.prefix}_*.sdf")
    print(f"Wrote energy/RMSD table to: {txt_path}")
    print("Index 0 is the minimum-energy conformer (within the selected set).")


if __name__ == "__main__":
    main()
