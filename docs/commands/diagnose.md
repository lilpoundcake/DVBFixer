# dvbfixer diagnose — structure quality report

[← command index](index.md) · [← README](../../README.md)

Inspects a PDB file and emits a plain-text per-residue findings report.
**Report-only** — never mutates the input. Users pick the appropriate
dvbfixer tool (`prepare` / `minimize` / `protonate` / `pull`) to
actually fix any issue.

Modelled after MolProbity's all-atom validation and BioLuminate's
Protein Report widget. Three families of checks run against every
input:

- **Structural integrity** — missing atoms, missing residues,
  missing terminals, coincident atoms, misplaced hydrogens, altLoc
  conflicts, chain breaks, insertion codes.
- **Chemistry / bond geometry** — valence violations, bond-length
  outliers, cis peptides, non-planar amides, Cα chirality.
- **Steric analysis** — all-atom clashes using OpenMM's neighbor
  search + van der Waals radii (or MolProbity's `probe` binary if
  it's on `PATH`).

## Usage

```bash
# Full report on stdout
dvbfixer diagnose input.pdb

# Write report to file
dvbfixer diagnose input.pdb -o report.txt

# Restrict to one category
dvbfixer diagnose input.pdb --only steric

# Show only ERROR findings (usable in CI / pre-commit gates)
dvbfixer diagnose input.pdb --severity ERROR

# Verbose (per-check timing + expanded atom detail)
dvbfixer diagnose input.pdb -v
```

## Exit codes

- **0** — no ERROR-level findings (structure is pipeline-safe).
- **1** — at least one ERROR-level finding remains after the severity
  filter. Suitable for shell gates:

  ```bash
  dvbfixer diagnose input.pdb --severity ERROR && \
      dvbfixer prepare input.pdb -o prepared.pdb
  ```

- **2** — argument or I/O error (input file not found, etc.).

## Severity levels

- **ERROR** — would break downstream tools if left as-is (missing
  atoms, coincident atoms, valence violations, hard clashes with
  ≥ 0.5 Å vdW overlap).
- **WARNING** — tolerated by the pipeline but suspect (moderate
  clashes 0.4 – 0.5 Å overlap, altLoc conflicts, chain breaks,
  near-cis peptides).
- **INFO** — noted for user review (insertion codes on antibody
  CDRs, cis-PRO — natural but worth flagging).

## Options

```
usage: dvbfixer diagnose [-h] [-o OUTPUT]
                         [--only {all,structural,chemistry,steric}]
                         [--severity {ERROR,WARNING,INFO}]
                         [--include-water]
                         [--format {text,json}] [-v]
                         input
```

See the auto-generated [reference page](../reference/diagnose.md) for
the full `--help` output.

Notable non-default behaviours:

- **Waters are excluded from chain-break AND steric checks by default.**
  Crystallographic waters generate massive noise (every ordered water
  triggers a chain break vs its neighbour). Pass `--include-water` to
  restore the old behaviour.
- **Multi-MODEL PDBs are analysed on MODEL 1 only.** A WARNING banner
  is emitted at the top of the report noting how many MODELs were
  detected. Use `dvbfixer split` if you need per-frame analysis.
- **JSON output** is available via `--format json` — a
  machine-readable list of findings suitable for CI gating or
  scripted post-processing.

## Sample output

```
================================================================================
dvbfixer diagnose — input.pdb
================================================================================
Loaded: 23 atoms, 2 residues, 1 chains

Structural integrity
----------------------------------------
  ERROR   B/SER126:HG            coincident with OXT (0.001 Å apart)
                                 fix: dvbfixer prepare (0.4.1 pre-strip re-places both)
  ERROR   B/SER126:HG            1.68 Å from parent OG (expected ~0.97 Å)
                                 fix: dvbfixer prepare (0.4.2 post-check re-places)

Chemistry / bond geometry
----------------------------------------
  ERROR   B/SER126:HG-OG         bond H-O: 1.677 Å (expected 0.960 Å, deviation 74.7%)
                                 fix: dvbfixer minimize

Steric analysis
----------------------------------------
  ERROR   B/SER126:OG            clashes with B/SER126:OXT — overlap 1.36 Å

Summary
----------------------------------------
  ERROR:   4 findings across 1 residues
  WARNING: 0 findings across 0 residues
  INFO:    0 findings

Suggested next step: `dvbfixer prepare input.pdb`
================================================================================
```

## Clash-detection engine

Two engines with the same output shape:

- **Python (default)** — `scipy.spatial.cKDTree.query_pairs` for the
  neighbor search, then a vdW-overlap test against the standard
  MolProbity radius table (H=1.20 Å, C=1.70 Å, N=1.55 Å, O=1.52 Å,
  S=1.80 Å, P=1.80 Å). Excludes 1-2, 1-3, and 1-4 bonded pairs so
  tight but chemically-valid rotamers aren't flagged. The bond graph
  is built from BOTH the OpenMM topology bond set AND distance-based
  inference (heavy-heavy ≤ 1.9 Å, X-H ≤ 1.3 Å, S-S ≤ 2.25 Å) — this
  catches HETATMs, glycans, and AMBER/CHARMM protonation variants
  whose real bonds OpenMM's PDBFile parser doesn't infer. Pairs that
  fit a hydrogen-bond envelope (polar H covalently bonded to N/O/S,
  paired with an N/O/S/F acceptor at 1.4 – 2.6 Å) are also skipped
  — MolProbity's `probe` classifies those as `hbond` contacts, not
  clashes.
- **MolProbity `probe`** — used automatically when the `probe` binary
  is on `PATH` (bundled with Phenix / MolProbity). Higher fidelity
  matching MolProbity's own reports. Falls back to the Python engine
  on any probe failure.

Both use the same clash thresholds — calibrated to MolProbity:

- **ERROR** — vdW overlap ≥ 0.5 Å.
- **WARNING** — vdW overlap 0.4 – 0.5 Å.
- Overlaps < 0.4 Å are below MolProbity's official clashscore floor
  and skipped.

## When to run this

- **Before every `dvbfixer prepare` / `minimize` / `zbs` run** on
  input you don't fully trust (RCSB downloads, upstream-tool output,
  hand-edited PDBs). Saves debugging a downstream crash.
- **As a pre-commit gate** for structural data checked into a repo:
  `dvbfixer diagnose file.pdb --severity ERROR` in the hook.
- **After suspicious pipeline output** to check whether a step
  produced something reasonable. Complementary to
  `dvbfixer.ffutils.geometry.repair_misplaced_hydrogens` which runs
  automatically inside `prepare` / `protonate` / `minimize`.

## See also

- [`prepare`](prepare.md) — auto-fixes coincident atoms + misplaced H
- [`minimize`](minimize.md) — energy-minimize to relax bond lengths /
  clashes
- [`pull`](pull.md) — targeted geometry fixes (SS bonds, glycosidic
  linkages)
- [BEST_PRACTICES.md](../../BEST_PRACTICES.md) — recipe for QA-first
  glycoprotein preparation
