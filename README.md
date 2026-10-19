<p align="center">
  <a href="https://github.com/lenhanpham/rmsd">
    <picture>
      <img src="resources/rmsd-logo.svg" alt="rmsd" style="width: 50%;">
    </picture>
  </a>
</p>
<p align="center">
  <a href="https://www.rust-lang.org/">
    <img src="https://img.shields.io/badge/Rust-1.85+-orange.svg" alt="Rust">
  </a>
  <a href="./LICENSE">
    <img src="https://img.shields.io/badge/License-MIT-green.svg" alt="License">
  </a>
  <a href="https://github.com/lenhanpham/rmsd/actions/workflows/build.yml">
    <img src="https://github.com/lenhanpham/rmsd/actions/workflows/build.yml/badge.svg" alt="Build Status">
  </a>
</p>

# rmsd 

`rmsd` is a program written in Rust that calculates how different two
geometrical structures are. It reads molecular geometries from XYZ or
PDB files and reports their **root-mean-square deviation (RMSD)** - a
single number, in the same length unit as the input coordinates
(typically Ångström), where `0` means identical geometries and larger
values mean larger average atomic displacements.

Given two structures A and B with N atoms at positions **p**ᵢ and
**q**ᵢ, the program evaluates

```text
RMSD = sqrt( Σᵢ |pᵢ − qᵢ|² / N )
```

after bringing both structures into the same frame: it centers them,
finds the optimal rotation, optionally reorders atoms that are listed
in different orders, and optionally searches reflections and axis
swaps. It works on single pairs (`rmsd a.xyz b.xyz`) as well as on
whole structure lists, producing a pairwise RMSD matrix for heatmap
visualisation (`rmsd --from-list list.txt`).

## Supported methods

Rotations: `kabsch` (default), `quaternion`, `none`.

Reorder methods: `none`, `hungarian`, `inertia-hungarian` (default),
`distance`, `brute` (reference only, exponential cost).

## About RMSD

As a rule of thumb, values below ~0.5 Å indicate nearly identical
structures, 1–2 Å a clearly related but relaxed geometry, and anything
much larger a different conformation or a bad atom correspondence.

Typical uses:

- **Conformer analysis** — quantify how different generated conformers
  of a molecule are (e.g. the 16 [5]dendralene geometries in
  `tests/5dendralene/`), cluster them, or pick diverse subsets.
- **Geometry optimization checks** — compare structures before and
  after relaxation, or from different optimizers.
- **Trajectory / ensemble analysis** — pairwise RMSD matrices feed
  directly into heatmaps and clustering of MD frames or docking poses.
- **Method benchmarking** — score predicted structures against
  reference geometries.

Because raw RMSD depends on the input frame and atom ordering, this
tool aligns structures first (see Theory): comparing without alignment
measures translation/rotation artifacts rather than chemistry.

## Theory

**Centering.** Both structures are translated to their centroid
(`Σ pᵢ / N`), or to the mass-weighted center of mass for
inertia-based methods.

**Rotation — Kabsch (default).** From the covariance `C = PᵀQ`, the
singular value decomposition `C = VSW` gives the optimal rotation
`U = VW`, with the last column of `V` flipped when
`det(V)·det(W) < 0` so the result stays right-handed
An optional weight vector produces the weighted variant, and
`kabsch_fit` additionally returns translated coordinates.

**Rotation — quaternion.** Horn's quaternion method
a 4×4 matrix accumulated from per-atom quaternion blocks is
diagonalized, and the eigenvector of the largest eigenvalue is the
optimal rotation. Numerically equivalent to Kabsch; offered as an
alternative rotation method.

**Atom reordering.** All reorderings work per chemical element, so
element multisets are always preserved:

| Method | Idea | When to use |
|---|---|---|
| `none` | atoms must already match in order | identical orderings |
| `hungarian` | optimal atom↔atom assignment minimizing summed distances (built-in Munkres implementation) | already-aligned frames, shuffled order |
| `inertia-hungarian` (default) | align inertia-tensor eigenvectors, then Hungarian over all 8 sign reflections, keeping the lowest-RMSD pairing | unaligned + shuffled (general case) |
| `distance` | match atoms sorted by distance from centroid | cheap heuristic |
| `brute` | exhaustive n! permutations per element | tiny systems only; reference |

**Reflections.** `--use-reflections` scans all 6 axis permutations × 8
sign reflections (48 trials) and keeps the best RMSD — this can match
mirror images but alters stereochemistry;
`--use-reflections-keep-stereo` skips enantiomeric trials.

## Install & build

Requires a recent stable Rust toolchain (edition 2024, rustc 1.85+):

```sh
git clone https://github.com/lenhanpham/rmsd.git
cargo build --release
./target/release/rmsd --version   # rmsd 0.1.0
```

Export the binary directory in your .bashrc if you want to use rmsd directly from your terminal 
```sh
export PATH=/path/to/rmsd/target/release:$PATH 
```

## Usage

```text
rmsd FILE_A FILE_B                 # pair mode: print one RMSD
rmsd -p FILE_A FILE_B              # pair mode: print rotated structure B
rmsd --from-list LIST              # matrix mode
rmsd LIST                          # matrix mode (single arg = list file)
rmsd --merge-csv OUT PART...       # merge partial matrix CSVs
```

Common flags: `-r/--rotation kabsch|quaternion|none`,
`-e/--reorder`, `--reorder-method …`, `--use-reflections`,
`--use-reflections-keep-stereo`, `--only-alpha-carbons` (PDB),
`-n/--ignore-hydrogen`, `--remove-idx/--add-idx`, `--format xyz|pdb`
(auto-detected, `.gz` transparent), `--threads`, `--max-mem`,
`--rows START:END`, `--out`, `--overwrite`, `--resume`. `--print` is pair-mode only.

## Examples

```sh
# Same water in a rotated/translated frame: ~0
$ rmsd tests/golden/water_a.xyz tests/golden/water_b.xyz
0.000000001682525212827008

# Two optimized ethane rotamers, atoms shuffled: reorder first
$ rmsd -e --reorder-method hungarian tests/golden/ethane_a.xyz tests/golden/ethane_b.xyz
0.45877449558576744

# Heavy atoms alone: the rotamers nearly coincide
$ rmsd -n -e --reorder-method hungarian tests/golden/ethane_a.xyz tests/golden/ethane_b.xyz
0.0009205572465018577

# Reflection search over two [5]dendralene conformers
$ rmsd --use-reflections tests/5dendralene/5dendralene_01.xyz tests/5dendralene/5dendralene_02.xyz
0.8926868464424043

# Rotated structure instead of a number
$ rmsd -p tests/golden/water_a.xyz tests/golden/water_b.xyz
3
Rotated 'tests/golden/water_b.xyz' to match 'tests/golden/water_a.xyz' with RMSD 0.00000000
O       0.00000000     -0.00000000      0.00000000
H       0.75700000      0.58600000      0.00000000
...
```

## Matrix mode

Compute every pair once and write a CSV heatmap table beside the list
file (`list-2.txt` → `list-2.csv`):

```sh
rmsd --from-list tests/5dendralene/list-1.txt
rmsd tests/5dendralene/list-2.txt   # single arg = list mode
```

List files accept one filename per line, whitespace-separated, or
comma-separated entries (all `tests/5dendralene/list-*.txt` styles).
Relative entries resolve against the list file's directory. The CSV has
a label header row + first column (names as listed), the lower triangle
filled including a `0.0` diagonal (6 decimals), and the upper triangle
blank — paste into Excel and apply a color-scale heatmap. All rotation /
reorder / reflection / filter flags apply to every pair.

Parallel + HPC runs — rows compute in parallel with Rayon (per-pair
results are bit-identical to serial runs):

```sh
# Headnode-safe defaults: half the cores, 4 GB molecule budget
rmsd --from-list big.txt
# Explicit resources (e.g. inside a queue allocation)
rmsd --from-list big.txt --threads 32 --max-mem 64GB
# Split rows across queue jobs, then merge (SLURM array example):
#   #SBATCH --array=0-3
#   rmsd --from-list big.txt --rows $((SLURM_ARRAY_TASK_ID*250)):$(...) --out part-$SLURM_ARRAY_TASK_ID.csv
rmsd --merge-csv full.csv part-0.csv part-1.csv part-2.csv part-3.csv
```

By default an existing output CSV is always archived beside the list file
as `<stem>_<YYYYMMDD_HHMMSS>[_<n>]<ext>` and the matrix is recalculated
from scratch. `--resume` instead validates the existing CSV (header and
row labels), discards any torn trailing line, and appends only the missing
rows; `--overwrite` rewrites from scratch without archiving.
Output streams in row blocks, so peak RAM stays at molecules + one
block, never the full matrix. Scale note: pairs grow as n(n+1)/2, so
million-structure sets must be chunked across many queue jobs —
one node cannot hold or compute the full matrix in reasonable time.

## Project layout

- `src/elements.rs` — Z/symbol/mass tables
- `src/geom.rs` — centroid, RMSD, Kabsch (+weighted), quaternion
- `src/inertia.rs` — center of mass, inertia tensor, principal axis
- `src/reorder.rs` — Hungarian (built-in Munkres), distance,
  inertia-hungarian, brute
- `src/reflect.rs` — 48 swap×reflection search
- `src/io/` — XYZ/PDB readers, XYZ writer, gzip
- `src/pipeline.rs` — filter → center → reorder/reflect → RMSD → transform
- `src/list.rs` — structure-list tokenizer + path resolution
- `src/matrix.rs` — pairwise lower-triangle compute + CSV serializer
- `src/bin/rmsd.rs` — `clap` CLI (pair + matrix modes)
- `tests/reference.rs` — 26 integration tests on real optimized molecules:
  16 [5]dendralene conformers (`tests/5dendralene/`) plus golden pairs
  (`tests/golden/`: ethane rotamers, rotated-frame C3H6, water, ALA PDB),
  all checked against reference values

## Verify

```sh
cargo test
cargo clippy --all-targets -- -D warnings
```

Reference values (same fixtures): dendralene Kabsch `1.49–1.99`,
ethane rotamers `0.947` (heavy atoms `9e-4`), rotated-frame C3H6 `~0`,
water `1.68e-09`, ALA PDB `0.01244` — all agree to ≤1e-9 absolute.
Known numerical sensitivity (documented in `tests/reference.rs`):
on symmetric tops (degenerate inertia moments, e.g. ethane/mirror),
`inertia-hungarian` depends on eigensolver details and `distance`
on tie-breaking, so those combos assert determinism/sanity rather than
exact golden values.
