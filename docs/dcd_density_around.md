# Description of `dcd_density_around.py`

## Purpose

Compute a **3D spatial density map** of molecule B (mol-B) around molecule A (mol-A) from a PBC-wrapped MD trajectory in DCD format. The output is an **OpenDX volumetric file** (viewable in VMD, PyMOL, etc.) representing how often mol-B beads appear in each voxel of a cube centered on mol-A.

This analysis answers **where mol-B tends to sit relative to mol-A.** This is done by repeatedly superimposing every mol-A copy in every frame onto a fixed reference orientation, applying the same rigid transform to nearby mol-B beads, and accumulating where those beads land. 

## Inputs

| Argument | Meaning |
|---|---|
| `pdb` | System PDB — defines the topology (all chains, all beads) |
| `dcd` | PBC-wrapped DCD trajectory |
| `native` | A single mol-A PDB used as the **superposition reference** |
| `--mol_a FIRST LAST` | Which chains (1-based range) are copies of mol-A |
| `--mol_b FIRST LAST` | Which chains are mol-B (default: everything not in mol-A) |
| `--range R` | Half-size of the analysis box in Å; defines the cubic region `[-R, R]³` |
| `--grid_size dx` | Voxel spacing in Å |
| `--output` | Output `.dx` file name |
| `--frames START END STRIDE` | Optional 1-based inclusive frame range with stride (default: all frames, stride 1) |
| `--movie` | Optional output PDB movie file |
| `--mol_a_fit_residues` | Optional comma-separated `res_seq` ranges (inclusive) selecting which mol-A residues are used for the superposition fit, e.g. `"2-10,15-20"`. Default: all mol-A atoms. |

----

## Setup (before the loop)

1. **System PDB parsing.** All chains are read. For each chain, two data structures are precomputed:
   - `slice` — flat bead-array index range for that chain, to extract its rows from DCD coordinate arrays.
   - `metadata` — per-bead `(name, res_name, chain_id, res_seq, ins_code)` tuples for PDB writing.

2. **Conservative radius for mol-B.** For each mol-B chain, `chain_max_radius` computes the maximum distance from the chain's geometric center to any of its beads (from static PDB coords). The **maximum** of these per-chain radii across all mol-B chains is used as a single global upper bound (`mol_b_global_radius`) for the spatial filter in Step 3, since mol-B may be flexible during simulation and per-chain static radii may underestimate the true extent.

3. **Native mol-A reference.** The native PDB is read and its bead coordinates are stored as an `(N_a, 3)` array, ready for `calc_rotation`. Without `--mol_a_fit_residues`, the native must have the same atom count as every mol-A copy in the system. With `--mol_a_fit_residues`, only the listed residues (matched by `res_seq`) are required in the native; any additional residues are silently ignored, so a small fragment PDB can serve as the alignment reference. In either case, the fit-subset atom count must match between native and every mol-A copy.

4. **Grid allocation.** A `(ngrid, ngrid, ngrid)` array of zeros is created where `ngrid = floor(2R / dx)`.


## Per-frame loop

For each sampled DCD frame in the requested sequence `START, START + STRIDE, START + 2*STRIDE, ...` up to `END`:

### For each mol-A copy (there can be multiple identical chains in the system)

**Step 1 — Unwrap mol-A.**
The raw PBC-wrapped coordinates of this mol-A copy are unwrapped using the **sequential MAXD algorithm**: beads are walked in order; if bead `i` jumps more than `MAXD = 25 Å` from bead `i-1`, the entire box vector is added or subtracted. This assumes connectivity is sequential and chain bonds never exceed 25 Å.

**Step 2 — Superimpose onto native.**
`calc_rotation(ref_np, mol_a_unwrp)` (from `fQCP`) computes the optimal rotation+translation that minimizes RMSD between the unwrapped mol-A copy and the native reference. It returns:
- `rmsd` — the fit quality
- `mat` — a **4×4 homogeneous transformation matrix** encoding the combined rotation and translation

This matrix is applied to mol-A itself (for the movie output) and to each nearby mol-B chain.

**Optional fit subset (`--mol_a_fit_residues`).** When supplied, only beads whose `res_seq` falls in one of the listed inclusive ranges are passed to `calc_rotation`. The resulting 4×4 matrix is still applied to the **full** mol-A copy and **full** mol-B chains — only the QCP fit itself uses the subset. This is useful when only part of mol-A is rigid enough to define a meaningful alignment (e.g., a structured core with flexible tails), or when only a fragment of mol-A has a reliable reference structure. In this mode, the **subset centroid** (rather than the full mol-A centroid) is used as the mol-A anchor in Step 3: under a subset fit, it is the subset — not the full mol-A — that maps onto its native-frame counterpart, and that is the sim-frame point that lands near the centre of the analysis box after the transform.

**Step 3 — Find relevant mol-B chains (conservative filter).**
For each mol-B chain:
- Its raw center is computed.
- **Minimum-image convention:** the shift that brings the mol-B center closest to the mol-A center under PBC is computed and applied to the entire chain.
- **Conservative skip condition:** the chain is skipped only if, in any dimension, the distance between centers exceeds `R + mol_b_global_radius` (the universal upper bound computed in Setup). That is: even if the chain were perfectly oriented toward the box, none of its beads could reach `[-R, R]³`. Chains that pass this filter are not guaranteed to have beads inside the box; those are filtered at binning time.

**Step 4 — Transform mol-B beads.**
The same `mat` (4×4 matrix from mol-A superposition) is applied to mol-B's shifted+unwrapped coordinates. This places mol-B in the **reference frame of the native mol-A** — the same frame for every snapshot, enabling meaningful accumulation.

**Step 5 — Bin beads.**
For each transformed mol-B bead, its position is converted to grid indices: `i = floor((x - origin) / dx)`. Only beads where all three indices fall within `[0, ngrid)` — i.e., truly inside `[-R, R]³` — increment the grid counter. This is the strict containment check.

**Optional movie output.** If `--movie` is given, one `MODEL` block is written per (frame, mol-A copy) pair, containing the transformed mol-A and all nearby mol-B chains.


## Normalization and output

After all sampled frames:
- `grid` holds raw bead counts (total hits per voxel summed across all frames and mol-A copies).
- It is divided by `total_samples × dx³` to give units of **beads/Å³** — a proper volumetric density.
- Written as an **OpenDX file**, which VMD or PyMOL can display as an isosurface or volumetric map overlaid on structure.
- The integral `density.sum() × dx³` is printed, which equals the **average number of mol-B beads inside the box per (frame, mol-A copy) sample**.
