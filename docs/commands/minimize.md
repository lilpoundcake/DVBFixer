# dvbfixer minimize — Energy Minimization with OpenMM

[← command index](index.md) · [← README](../../README.md)

Energy-minimizes a PDB structure with OpenMM using selective restraints. **By default minimizes the whole system** (protein + sugars + ligands) using the resolved `--ff` (short-name aliases + auto-detection — see [force-fields.md](../force-fields.md)). When GLYCAM residues are present, the FF is upgraded to AMBER14 + GLYCAM_06j-1 automatically; `add_glycam_bonds(positions=...)` populates intra-residue bonds for GLYCAM templates + protein-glycan peptide bonds + distance-based sugar-sugar glycosidic bonds, and `ignoreExternalBonds=True` is used for N-linked glycan junctions. For arbitrary unknown ligands (drug molecules, cofactors) that lack a template in the resolved FF, add `--parametrize-ligands` and dvbfixer runs GAFF2 + AM1-BCC via antechamber and registers a real OpenMM template per ligand (cached under `~/.cache/dvbfixer/lig_params/`). Pass `--strip-heterogens` for the strip-and-splice flow (protein-only minimization with HETATM coords restored verbatim). NLN/OLS/OLT names from the input PDB are snapshotted and restored just before the final write — defensive belt-and-braces for the strip-and-splice fallback. Reads a `.dat` file (from `dvbfixer model` + `dvbfixer prepare`) to apply different restraint strengths to original vs newly added atoms; heterogens not in `.dat` get no restraint and relax freely. All HIS residues are automatically renamed to explicit variants (HIE/HID/HIP). By default keeps existing hydrogens from input; with `--rebuild-h`, strips and re-adds via OpenMM. Detects AMBER protonation names (HIE/GLH/CYX etc.) from the raw PDB and passes them as `variants` to `addHydrogens`. Calls `fix_atom_hetatm_records` on the output so protonation variants and NLN/OLS/OLT are written as ATOM records.

For BioLuminate-style refinement of arbitrary ligands/sugars **without per-ligand parametrization**, pass `--xtb-refine` (xtb GFN-FF universal force field) or `--obminimize-refine` (OpenBabel MMFF94/UFF/GAFF) — both auto-type any organic molecule from connectivity rules. Combine with `--refine-heterogens-only` to keep the AMBER-quality protein frozen and refine only the glycan/ligand geometry. Anchor residues (e.g. ASN of an N-linked glycan) are included as frozen atoms so protein-glycan bond geometry is preserved.

> **Note**: `--obminimize-refine` gives sharper glycosidic bond geometry (e.g. ASN-NAG ~1.48 Å). `--xtb-refine` is slower and the v6.7.1 build in conda-forge has a known `$fix` bug that can stretch the linkage to ~1.66 Å. Prefer obminimize until xtb 6.8+ ships.

## Three-Tier Restraint System

| Atom category | Force constant | Purpose |
|---------------|----------------|---------|
| Original heavy atoms | 100 kcal/mol/A^2 | Keep resolved structure in place |
| New backbone (N, CA, C, O, CB) | 5 kcal/mol/A^2 | Maintain reasonable loop geometry |
| New sidechain + all hydrogens | 0 (free) | Full relaxation |

Minimization runs in two phases:
1. Full restraints (1000 iterations)
2. Restraints reduced 10x (1000 iterations)

## Usage

```bash
# Minimize with .dat restraint info
dvbfixer minimize input_prepared.pdb --dat input_prepared.dat -v
# Outputs: input_prepared_minimized.pdb

# Vacuum minimization (no solvent box)
dvbfixer minimize input.pdb --no-solvent

# Without .dat — all atoms get strong restraints
dvbfixer minimize input.pdb

# BioLuminate-style: AMBER protein + xtb GFN-FF refinement for any heterogen
# (no per-ligand parametrization, works for arbitrary sugars/ligands)
dvbfixer minimize glycoprotein.pdb --xtb-refine --refine-heterogens-only

# Faster alternative using OpenBabel UFF (seconds vs minutes for xtb)
dvbfixer minimize glycoprotein.pdb --obminimize-refine --refine-heterogens-only
```

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o`, `--output` | `<input>_minimized.pdb` | Output minimized PDB |
| `--dat` | `<input>.dat` | Restraint data file from `dvbfixer prepare` |
| `--ph` | 7.0 | pH for hydrogen addition if needed |
| `--ff` | `auto` | Force field. Accepts a short name (`auto`, `amber`, `amber+glycam`, `charmm`, …) or explicit OpenMM XML paths. See [force-fields.md](../force-fields.md). |
| `--padding` | 1.0 | Solvent padding in nm |
| `--restraint-k` | 100.0 | Strong restraint constant (kcal/mol/A^2) |
| `--weak-k` | 5.0 | Weak restraint constant for new backbone (kcal/mol/A^2) |
| `--max-iter` | 1000 | Max minimization iterations per phase |
| `--rebuild-h` | off | Strip and re-add hydrogens via OpenMM (default: keep existing) |
| `--strip-heterogens` | off (default: keep) | Strip heterogens before parametrization, splice coords back — protein-only mode. **Warning**: heterogens are restored at their INPUT coords while the protein has moved during minimization; the protein-ligand interface (H-bonds, contacts) may end up strained. For proper interface geometry use `--parametrize-ligands` instead. |
| `--no-solvent` | off | Minimize in vacuum |
| `--xtb-refine` | off | Post-pass: refine geometry with xtb GFN-FF universal force field (auto-parametrizes any organic molecule, no templates needed) |
| `--xtb-cycles` | 200 | Max xtb optimization cycles |
| `--obminimize-refine` | off | Post-pass: refine geometry with OpenBabel obminimize (faster than xtb, UFF / MMFF94 / GAFF) |
| `--obminimize-ff` | UFF | OpenBabel force field. UFF is default — handles N-glycosidic linkages correctly. MMFF94s mistypes anomeric C as sp2 (gives ~120° instead of ~109° angles around the C1-N bond) |
| `--obminimize-steps` | 500 | OpenBabel minimization steps |
| `--refine-heterogens-only` | off | With `--xtb-refine`/`--obminimize-refine`: refine only heterogen residues (protein frozen). BioLuminate-style ligand-only minimization. **Caveat**: only the ligand's INTERNAL geometry gets refined — the protein-ligand INTERFACE is NOT relaxed, so any pre-existing clash there persists. Drop the flag for whole-system refinement when the interface matters. |
| `--platform` | auto | OpenMM platform (CPU, CUDA, OpenCL, Reference) |
| `--rename` | off | Rename non-canonical residues (AMBER/CHARMM) to standard names before processing |
| `--parametrize-ligands` | off | For each heterogen residue with no template in the resolved `--ff`, run GAFF2 + AM1-BCC via antechamber and register the result as an OpenMM template. Cached under `~/.cache/dvbfixer/lig_params/` (override with `$DVBFIXER_LIG_CACHE`). Requires AmberTools (`antechamber`/`parmchk2`) and `openmmforcefields`. See [force-fields.md](../force-fields.md). |
| `-v`, `--verbose` | off | Print detailed progress |

## Recommended Workflow

```bash
# 1. Fix structure
dvbfixer prepare input.pdb -v
# -> input_prepared.pdb + input_prepared.dat

# 2. Inspect in PyMOL/VMD/ChimeraX, edit .dat if needed

# 3. Minimize
dvbfixer minimize input_prepared.pdb --dat input_prepared.dat -v
# -> input_prepared_minimized.pdb
```

## See also

- [`prepare`](prepare.md) — produces the `.dat` file `minimize` consumes
- [`protonate`](protonate.md) — assign protonation state names between two `minimize` passes
- [BEST_PRACTICES.md](../../BEST_PRACTICES.md) — xtb / obminimize refinement notes
