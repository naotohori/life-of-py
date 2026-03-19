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
| `--frames START END` | Optional frame range (default: all) |
| `--movie` | Optional output PDB movie file |

----

## Setup (before the loop)

1. **System PDB parsing.** All chains are read. For each chain, two data structures are precomputed:
   - `slice` — flat bead-array index range for that chain, to extract its rows from DCD coordinate arrays.
   - `metadata` — per-bead `(name, res_name, chain_id, res_seq, ins_code)` tuples for PDB writing.

2. **Conservative radius per mol-B chain.** For each mol-B chain, `chain_max_radius` computes the maximum distance from the chain's geometric center to any of its beads (from static PDB coords). This is used later as a spatial filter.

3. **Native mol-A reference.** The native PDB is read, its bead coordinates stored as `ref_mol_a_F` in Fortran-order `(3, N_a)` shape, ready for `calcrotation`.

4. **Grid allocation.** A `(ngrid, ngrid, ngrid)` array of zeros is created where `ngrid = floor(2R / dx)`.


## Per-frame loop

For each DCD frame within the requested range:

### For each mol-A copy (there can be multiple identical chains in the system)

**Step 1 — Unwrap mol-A.**
The raw PBC-wrapped coordinates of this mol-A copy are unwrapped using the **sequential MAXD algorithm**: beads are walked in order; if bead `i` jumps more than `MAXD = 20 Å` from bead `i-1`, the entire box vector is added or subtracted. This assumes connectivity is sequential and chain bonds never exceed 20 Å.

**Step 2 — Superimpose onto native.**
`calcrotation(ref_mol_a_F, mol_a_F)` (from an external `fQCP` library) computes the optimal rotation+translation that minimizes RMSD between the unwrapped mol-A copy and the native reference. It returns:
- `rmsd` — the fit quality
- `mat` — a **4×4 homogeneous transformation matrix** encoding the combined rotation and translation

This matrix is applied to mol-A itself (for the movie output) and to each nearby mol-B chain.

**Step 3 — Find relevant mol-B chains (conservative filter).**
For each mol-B chain:
- Its raw center is computed.
- **Minimum-image convention:** the shift that brings the mol-B center closest to the mol-A center under PBC is computed and applied to the entire chain.
- **Conservative skip condition:** the chain is skipped only if, in any dimension, the distance between centers exceeds `R + chain_max_radius`. That is: even if the chain were perfectly oriented toward the box, none of its beads could reach `[-R, R]³`. Chains that pass this filter are not guaranteed to have beads inside the box; those are filtered at binning time.

**Step 4 — Transform mol-B beads.**
The same `mat` (4×4 matrix from mol-A superposition) is applied to mol-B's shifted+unwrapped coordinates. This places mol-B in the **reference frame of the native mol-A** — the same frame for every snapshot, enabling meaningful accumulation.

**Step 5 — Bin beads.**
For each transformed mol-B bead, its position is converted to grid indices: `i = floor((x - origin) / dx)`. Only beads where all three indices fall within `[0, ngrid)` — i.e., truly inside `[-R, R]³` — increment the grid counter. This is the strict containment check.

**Optional movie output.** If `--movie` is given, one `MODEL` block is written per (frame, mol-A copy) pair, containing the transformed mol-A and all nearby mol-B chains.


## Normalization and output

After all frames:
- `grid` holds raw bead counts (total hits per voxel summed across all frames and mol-A copies).
- It is divided by `total_samples × dx³` to give units of **beads/Å³** — a proper volumetric density.
- Written as an **OpenDX file**, which VMD or PyMOL can display as an isosurface or volumetric map overlaid on structure.
- The integral `density.sum() × dx³` is printed, which equals the **average number of mol-B beads inside the box per (frame, mol-A copy) sample**.
