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

## How it works
Energy-minimizes with OpenMM using selective restraints from the `.dat` file. Three tiers: original heavy=strong (100), new backbone=weak (5), new sidechain+H=free. Two phases: full restraints then 10x reduced. **Default: keeps heterogens and minimizes the whole system** (protein + sugars + ligands). The `--ff` is resolved through `ffutils.resolve_ff` (short-name aliases + auto-detection — see the "Shared FF selection" section below); when GLYCAM residues are detected the alias upgrades to `amber+glycam` automatically. `create_forcefield_with_openff` loads the resolved XMLs and prunes GLYCAM sugar/NA templates that would fuzzy-match PDB-named sugars to the wrong entry. Then `acpype_export.add_glycam_bonds(positions=...)` populates intra-residue bonds + protein-glycan peptide bonds + sugar-sugar glycosidic bonds (distance-based), and `ignoreExternalBonds=True` is used for N-linked glycan junctions. Heterogens not in `.dat` get no restraint → free to relax (matches `transplant.py --relax`). Pre-solvent createSystem check passes `residueTemplates` (built from CYX/HIE/HID/HIP/LYN auto-detection) to avoid CYM/CYX false positives. `main()` snapshots NLN/OLS/OLT names from the raw input PDB at startup and restores them on the final topology just before `PDBFile.writeFile`. Calls `fix_atom_hetatm_records` on the output.

**`--parametrize-ligands`** (opt-in): for each heterogen residue that has no template in the resolved base FF (and isn't standard protein/water/ion/GLYCAM/CHARMM), extracts the residue's atoms + coords via OpenBabel → SDF, wraps in `openff.toolkit.Molecule`, and hands the list to `openmmforcefields.generators.GAFFTemplateGenerator` (real AMBER GAFF2 + AM1-BCC via antechamber under the hood). The generator is passed to `create_forcefield_with_openff` via `extra_generators=` and registered on the ForceField before `createSystem`. Cached on disk under `~/.cache/dvbfixer/lig_params/gaff_ligands.json` (override with `$DVBFIXER_LIG_CACHE`). Implementation in `src/dvbfixer/lig_params.py`. Same limitation as any per-molecule ligand FF: cross-residue bonds between two ligand residues get no parameters (use `dvbfixer convert` to bring glycans into a template-supported scheme). SMIRNOFF was previously used in this slot but never actually parameterised cross-residue glycan bonds — removed in favour of GAFF2 via GAFFTemplateGenerator. `_extract_residue_sdf` (since 0.7.9) builds a heavy-atom-only sub-molecule first (`ConnectTheDots()` + `PerceiveBondOrders()`), then re-attaches hydrogens at their existing positions with forced single bonds — running OpenBabel's whole-molecule bond perception directly on an already-hydrogenated residue badly miscalled bond orders. `dvbfixer.ffutils.ligand_valence` supplies known-alkene overrides (e.g. DAN's ring C2=C3) and a connectivity-based ionizable-group detector (carboxylate/sulfonate/phosphate) so the resulting GAFF2 template gets the correct bond graph and net charge instead of an unfillable radical.

`--strip-heterogens` opt-in mode: strip-and-splice path — heterogens removed before parametrization, minimized protein-only with original `--ff`, HETATM coords restored by `(chain, resid, atomname)` matching.

**Strip-and-splice fallback path** (automatic when full-system GLYCAM parametrization fails on PDB-named sugars like NAG/BMA/MAN without GLYCAM templates): `_rigid_track_glycan_trees` runs after the protein-only minimize, doing Kabsch tracking of each glycan tree to follow the post-minimize protein anchor. Computes ideal trans-amide C1 position via `u_C1 = cos(120°)*u_CG + sin(120°)*(-u_perp_toward_OD1)` (places C1 in the amide plane on the OPPOSITE side of OD1 across the CG-ND2 axis — canonical E,Z amide). Also snaps HD21 to the canonical cis-OD1 position. Heterogen subsystem extraction (`_extract_heterogen_subsystem`) freezes only the protein BACKBONE of the anchor residue (`N/CA/HA/CB/HB2/HB3/C/O/H` for ASN; `N/CA/HA/CB/HB2/HB3/CG2/HG21/HG22/HG23/C/O` for THR) — the amide group (CG, OD1, ND2, HD21) and the heterogen linkage atom (NAG C1) are LEFT FREE so the subsequent xtb/obminimize refinement pass actually minimizes the ND2-C1 bond, angle, and amide planarity. **Critical**: neither the pre-solvent fallback nor the post-createSystem fallback calls `_restore_glycosylated_h` — it would overwrite the post-min amide (CG/OD1/ND2/HD21) with stale prep coords, breaking the bond to the rigid-tracked glycan. The rigid-tracker owns the geometry.

**BioLuminate-style refinement passes** (for arbitrary heterogens, no per-ligand parametrization):
- `--xtb-refine` runs xtb GFN-FF (universal force field, ~minutes for 400-atom glycan tree) via subprocess on a temp XYZ file. `_run_xtb()` writes XYZ, runs `xtb input.xyz --opt --gfnff --cycles N --norestart`, parses `xtbopt.xyz`.
- `--obminimize-refine` runs OpenBabel obminimize (MMFF94/UFF/GAFF, seconds). `_run_obminimize()` writes PDB, runs `obminimize -ff <FF> -n <steps> in.pdb`, parses stdout PDB. UFF works for any atom; MMFF94 may fail with "USING EMPIRICAL RULE / could not find van der Waals parameters" on atypical sugar atoms.
- `--refine-heterogens-only` extracts a sub-topology via `_extract_heterogen_subsystem()` that includes: (a) all heterogen residues + their internal bonds, (b) the FULL anchor residues (e.g. ASN that bonds to NAG via ND2 — entire ASN gets included, not just ND2, otherwise xtb perceives a lone N and pulls the C-N bond to ~1 Å), (c) the cross-residue bond preserved. After refinement, anchor atoms are frozen via `$fix atoms:` xcontrol block (xtb) or `OBFFConstraints.AddAtomConstraint()` (obminimize); only heterogen heavy atoms move. Refined coords spliced back into full topology; anchor coords kept verbatim (don't overwrite OpenMM-minimized protein).
- **xtb v6.7.1 `$fix` bug**: ~0.1 Å drift on nominally-frozen atoms. Acceptable but the glycosidic bond may stretch slightly. Fixed in v6.8 (not yet in conda-forge). obminimize via `OBFFConstraints` is rock-solid by comparison.
- **obminimize Python API path** (`_run_obminimize_pybel`): when frozen atoms are needed, switches from CLI to pybel. `obff.Setup(mol.OBMol, constraints)` + `obff.SteepestDescent(steps)` + `obff.GetCoordinates(mol.OBMol)`. CLI obminimize has NO freeze flag.
- **Auto-switch for big systems**: when the user passes `--xtb-refine` or `--obminimize-refine` without `--refine-heterogens-only`, and topology has >5000 atoms, auto-switches to heterogens-only with an INFO note. Avoids the OOM (UFF) / hours-long (xtb) trap on whole-protein invocations.
- **obminimize MMFF94 → UFF auto-fallback**: MMFF94/MMFF94s often lacks VdW params for atypical sugar atoms ("COULD NOT FIND VAN DER WAALS PARAMETERS"). `_run_obminimize` detects this in stderr/stdout and silently retries with UFF (covers any element).
- **`--obminimize-ff` default is UFF** (not MMFF94s). Rationale: MMFF94/MMFF94s mistypes the anomeric C of N-glycosidic linkages (ASN ND2 - NAG C1) as sp2. Result: angles around C1 settle at ~120° instead of the correct ~109° tetrahedral. UFF uses geometry-from-atomic-number typing and gives proper sp3 geometry (ND2-C1-C2 ≈ 110°, ND2-C1-O5 ≈ 108°).
- **Frozen anchors include the heterogen-side linkage atom** (e.g. NAG C1) in addition to the full protein anchor residue (full ASN). UFF/MMFF lack the amide-planarity term, so without freezing C1 the linkage carbon can rotate ~15° out of the amide plane during refinement. Freezing it preserves the OpenMM-AMBER interface geometry; the rest of the sugar refines freely.
- **`_restore_glycosylated_h` post-AMBER**: legacy strip-and-splice fallback renames NLN→ASN before AMBER parametrization. AMBER doesn't know about the glycan and may flip HD21/HD22 to the wrong side of ND2 (collision with linked C1, broken amide planarity). After minimize+splice, restore only the **amide group atoms** (CG/OD1/ND2/HD21 for NLN; OG for OLS; OG1 for OLT) from prep positions. CB and HB2/HB3 are NOT restored — they get proper sp3 tetrahedral geometry from the AMBER minimization (HB2-CB-HB3 ≈ 107.5°, HB*-CB-CA/CG ≈ 108-110°). Restoring more than the amide breaks methylene angles; restoring less breaks amide planarity.
- `_find_binary()` looks for `xtb`/`obminimize` first in PATH, then in `os.path.dirname(sys.executable)` — handles direct-executable env invocation where the env's bin dir isn't on PATH.
- Both engines run AFTER OpenMM minimize, after solvent stripping, before final PDB write.

**Known bug fixed**: RDKit's `AddHs(addCoords=True)` sometimes leaves an H atom at origin (0,0,0) when 3D placement fails. `prepare.add_heterogen_h_via_rdkit` now detects this and offsets the H from its parent atom by ~1 Å — otherwise xtb crashes on "very short distance" between zero-coord atoms. All HIS residues renamed to explicit variants (HIE/HID/HIP) before any OpenMM operation. Default: keeps existing hydrogens from input; if any residues are missing H (e.g. from mutation), uses PDBFixer to add them (OpenMM's `addHydrogens` can't handle residues with no H). `--rebuild-h` strips existing H, runs PDBFixer to fix missing heavy atoms and terminal atoms (OXT), then re-adds correct H via OpenMM. Detects AMBER protonation names (HIE/GLH/CYX etc.) from raw PDB text before OpenMM normalizes them, and passes them as `variants` to `addHydrogens` so correct protonation hydrogens are added (e.g. HE2 for GLH). Also reads `variant_overrides` from `.dat` file (saved by prepare) to recover protonation info even when output PDB has standard names.

**Important:** OpenMM's `PDBFile` reader normalizes AMBER names (GLH→GLU, HIE→HIS, CYX→CYS). The raw PDB must be read first with `_read_amber_renames()` to capture original names before loading with PDBFile.

**Known issue:** `.dat` stores chain IDs from PDBFixer. External tools may reassign chain IDs, breaking `.dat` matching.
