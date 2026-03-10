# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

dvbfixer is a Python package providing CLI tools for preparing PDB (Protein Data Bank) structural biology files. Installed as a single `dvbfixer` command with subcommands: `split`, `renumber`, `model`, `pull`, `prepare`, `minimize`, `protonate`, `rename`, `top`, `transplant`, `puppet`, `zbs`.

## Package Structure

```
src/dvbfixer/
  cli.py          — unified CLI entry point (dvbfixer <command>)
  split_chains.py — empirical chain splitting
  renumber.py     — SEQRES-based residue renumbering
  model.py        — Modeller-based loop/gap rebuilding
  pull.py         — bond pulling via OpenMM partial minimization
  prepare.py      — PDBFixer: add missing atoms/residues
  minimize.py     — OpenMM energy minimization with selective restraints
  protonate.py    — PROPKA3 pKa-based protonation
  rename.py       — rename non-canonical residues to standard names
  top.py          — GROMACS topology (.itp/.top) generation
  rtp_parser.py   — GROMACS force field file parsers (RTP/ARN/R2B/TDB/ATP)
  transplant.py   — transplant molecules between PDB structures
  puppet.py       — strip PDB to backbone-only polyglycine
  zbs.py          — full pipeline (renumber->model->prepare->minimize->protonate->minimize)
  acpype_export.py — ACPYPE-based GROMACS topology export (AMBER+GLYCAM)
  ffutils.py      — shared force field utilities (residue sets, OpenFF setup)
  pdbutils.py     — shared PDB utilities (CONECT remapping, serial maps)
FF/               — bundled GROMACS force field directories (AMBER, CHARMM36)
```

## Environment & Installation

```bash
micromamba create -f environment.yml
micromamba run -n dvbfixer pip install -e .
```

After installation, `dvbfixer` is available as a CLI command:
```bash
micromamba run -n dvbfixer dvbfixer <command> [args]
```

Key packages: numpy, OpenMM 8.4, PDBFixer 1.12, scipy, PROPKA 3.5, Modeller 10.8, BioPython.

**Modeller requires a license key** from https://salilab.org/modeller/registration.html (free for academics). Set in `<env>/lib/modeller-10.8/modlib/modeller/config.py`.

## Subcommands

### dvbfixer split
Splits chains in PDB files lacking chain IDs (e.g. GROMACS output). Water and ions are stripped before chain detection to prevent false breaks (`--keep-water` re-appends them). Three detection criteria:
1. Residue number backward jump (insertion codes handled — equal resSeq with different iCode is NOT a break)
2. C->N peptide bond distance > 2.5 A (protein residues only, from `STANDARD_RESIDUES` set)
3. Nearest-atom gap > 15 A (fallback for non-protein: sugars, ligands)

### dvbfixer renumber
Renumbers residues by aligning ATOM records to SEQRES via subsequence matching. Removes insertion codes (e.g. Kabat 100A-J -> sequential). Updates **all** PDB sections: ATOM, HETATM, TER, HELIX, SHEET, SSBOND, LINK, CISPEP, HET, DBREF, SEQADV, CONECT, REMARK 465/500/610. Each section has specific column positions — see the `update_*` functions and `remap_resid()` helper.

### dvbfixer model
Rebuilds missing loops/gaps using Modeller's LoopModel. Takes SEQRES (or --fasta) as complete sequence, aligns to ATOM records via `align2d`, runs loop modeling with configurable MD refinement. Non-protein chains (glycans, ligands) are included in the Modeller pipeline via `env.io.hetatm=True` with `'.'` (BLK) entries in the target PIR sequence so Modeller keeps them through loop modeling. Post-processing restores: (1) original chain IDs from Modeller's A,B,C,..., (2) original residue numbering with insertion codes using the alignment (template positions get original (resSeq, iCode), gaps get interpolated numbers). `--no-terminal` trims N/C terminal missing residues from the target sequence. Terminal alignment is auto-fixed after `align2d` to prevent misplaced terminal gaps. Writes a `.dat` file recording all atoms in rebuilt gap residues — prepare merges this with its own additions. Water removed by default (`--keep-water` to preserve).

### dvbfixer pull
Pulls atoms together to form bonds (disulfide bridges, glycosidic bonds) using OpenMM partial minimization. Supports multiple `--bond` specifications. Protein-protein bonds use `CustomBondForce`; protein-HETATM bonds use `CustomExternalForce` toward the fixed HETATM position. Atoms within `--radius` of bond endpoints are free to move (mass=0 freezing for the rest). Auto-removes conflicting hydrogens (CYS HG for disulfides, ASN HD22 for glycosidic bonds). Pre-pull validation checks valence (bond count vs `MAX_BONDS`) and bond type reasonableness for the pulling residues. Post-pull validation checks convergence (distance vs target), bond length range, and steric clashes within pulling residues.

### dvbfixer prepare
Runs PDBFixer to add missing residues, heavy atoms, and hydrogens. Strips H before PDBFixer processing to ensure clean template matching (H re-added by `addMissingHydrogens`). Writes a `.dat` file (JSON) recording which atoms were added. Merges upstream `.dat` (from model step) if present — atoms from gap rebuilding carry through to minimize. Detects glycosylated residues from CONECT records and removes extra hydrogens (e.g. ASN HD22 on glycosylated ND2). Supports `--mutate CHAIN:RESNUM:NEW_AA` for point mutations (can be used multiple times). `--rename` canonicalizes non-standard residue names before processing.

### dvbfixer minimize
Energy-minimizes with OpenMM using selective restraints from the `.dat` file. Three tiers: original heavy=strong (100), new backbone=weak (5), new sidechain+H=free. Two phases: full restraints then 10x reduced. Non-protein residues (glycans, ligands) are stripped before parametrization and restored with original coordinates afterward. Default: keeps existing hydrogens from input (adds via OpenMM only if none present). `--rebuild-h` strips existing H, runs PDBFixer to fix missing heavy atoms and terminal atoms (OXT), then re-adds correct H via OpenMM. Detects AMBER protonation names (HIE/GLH/CYX etc.) from raw PDB text before OpenMM normalizes them, and passes them as `variants` to `addHydrogens` so correct protonation hydrogens are added (e.g. HE2 for GLH).

**Important:** OpenMM's `PDBFile` reader normalizes AMBER names (GLH→GLU, HIE→HIS, CYX→CYS). The raw PDB must be read first with `_read_amber_renames()` to capture original names before loading with PDBFile.

**Known issue:** `.dat` stores chain IDs from PDBFixer. External tools may reassign chain IDs, breaking `.dat` matching.

### dvbfixer protonate
Runs PROPKA3 for pKa prediction, renames residues to AMBER protonation names (HID/HIE/HIP, ASH, GLH, CYM, CYX, LYN) based on target pH. `--no-hydrogens` for text-based rename only (used in pipeline). Default mode strips H and re-adds via OpenMM with `variants` parameter.

### dvbfixer rename
Text-based rename of non-canonical residue names to standard PDB names. Converts AMBER protonation (HIE/HID/HIP→HIS, ASH→ASP, GLH→GLU, CYX/CYM→CYS, LYN→LYS), CHARMM (HSD/HSE/HSP→HIS), and MSE→MET. Also available as `--rename` flag on `renumber`, `prepare`, `minimize`, and `pull` tools. Uses `CANONICAL_MAP` dict in `rename.py`.

### dvbfixer top
Generates self-contained GROMACS `.top` topology files directly from PDB by parsing force field RTP/ARN/R2B/TDB files in Python (no `pdb2gmx`). Supports AMBER99SB-ILDN and CHARMM36 force fields (bundled in `FF/` directory). The output `.top` file contains all FF parameters inlined (`[ defaults ]`, `[ atomtypes ]`, `[ bondtypes ]`, bonded params, water moleculetype, ion moleculetypes) — no external FF directory needed in the GROMACS working directory.

**Topology generation algorithm:**
1. Resolve PDB residue names to RTP names via `CANONICAL_MAP` + R2B mapping
2. Map PDB atom names to FF atom names via ARN file
3. Build atom list from RTP `[ atoms ]` with types, charges, masses
4. Build bond graph from RTP `[ bonds ]`, resolving inter-residue refs (`-C`, `+N`)
5. Enumerate angles (all 3-atom paths), dihedrals (all 4-atom paths), 1-4 pairs
6. Copy impropers and CMAP (CHARMM) from RTP
7. Apply terminal patches (CHARMM: NH3+/COO- from TDB files; AMBER: separate NXXX/CXXX RTP entries)

**Glycan support (CHARMM):** Parses `carb.rtp`, maps PDB sugar names (NAG, BMA, MAN, GAL, FUC, SIA) to CHARMM RTP names via `PDB_TO_CARB`. Detects glycosidic bonds from C1-O distances < 2.0 A. Builds glycan trees via BFS. Handles charge redistribution (HO removal → charge to O) and atom type changes (OC311→OC3C61 at linkage sites). Inter-residue glycan angles use ftype=1 (harmonic, no Urey-Bradley).

**Inter-chain SS bonds:** Written to separate `interchain_ss.itp` with `[ intermolecular_interactions ]`. Must be appended to `.top` after `gmx solvate`/`genion` (cannot appear before `[ molecules ]`).

**Output:** Self-contained `topol.top` (all FF params + chain moleculetypes + water + ions inlined), `posre_*.itp` position restraint files (separate, used with `#ifdef POSRES`), and `conf.pdb` with FF atom names. No separate chain `.itp` files — everything is in the `.top`. FF directory is only used at build time for RTP parsing, not needed at runtime.

**Key writers:** `_write_moleculetype(f, chain_top, bonded_types)` writes a chain moleculetype section to an open file handle. `write_top()` inlines FF content via `_read_ff_content()` (strips preprocessor directives), `_parse_defaults()` (extracts `[ defaults ]`), and `_write_water_topology()` (extracts rigid settles version from water .itp).

**Key data structures:** `ChainTopology` dataclass holds atoms, bonds, pairs, angles, dihedrals, impropers, cmap per chain. `TopologyBuilder` class loads FF data, builds chains, writes output. `rtp_parser.py` provides `parse_rtp()`, `parse_r2b()`, `parse_arn()`, `parse_tdb()`, `parse_atomtypes()`.

**Atom name matching (`_match_atom_names`):** Multi-pass algorithm: (0) exact match, (0a) ARN reverse mapping (e.g. HN→H for CHARMM), (1) strip trailing digits, (2) common H renames, (3) numbered variants, (4) singleton numbered atoms (HG1→HG). Validates all RTP heavy atoms are present; warns on extra/missing atoms.

**`--acpype` mode:** Alternative to RTP-based topology. Uses `acpype_export.py` shared module: OpenMM (AMBER14+GLYCAM) → ParmEd → ACPYPE → GROMACS `.top`/`.gro` with `[ pairs_nb ]` for mixed 1-4 scaling. Ignores `--ff`/`--water`/`--ignh`/`--merge` flags. Respects `--ss` for explicit disulfide bonds.

### dvbfixer transplant
Transplants molecules from a graft PDB into an acceptor PDB, with optional AMBER+GLYCAM energy minimization. Designed for the GLYCAM glycoprotein workflow: extract glycosylation site residues → submit to GLYCAM-Web → transplant back with glycans attached.

**Graft workflow** (`--donor` + `--graft`):
1. `--donor`: original protein residues extracted from acceptor (e.g. ASN307/SER308/THR309). Identifies which residues to replace and provides CA atoms for alignment.
2. `--graft`: GLYCAM output with renamed protein residues (NLN/OLS/OLT) + glycan trees (UYB/4YB/VMB etc.)
3. Kabsch superposition aligns donor→acceptor; same transform applied to graft
4. Donor residues removed from acceptor, graft protein residues inserted at correct positions, glycans appended
5. CONECT records remapped via (chain, resseq, atomname) identity matching (not serial numbers, which collide between sources)

**Simple transplant** (`--donor` + `--select`, no `--graft`): copies selected chains/residues from donor to acceptor with alignment. Auto-remaps chain IDs on collision.

**`--relax` minimization** with AMBER14 + GLYCAM_06j-1 force fields:
- Preprocessing: CYS→CYX for disulfide bonds (detected from CONECT), strip H from GLYCAM protein residues (NLN/OLS/OLT), remove spurious terminal atoms (OXT/H2/H3) from mid-chain graft residues
- Adds intra-residue bonds for GLYCAM residues from FF templates (OpenMM PDBFile doesn't infer bonds for non-standard residues)
- Adds peptide bonds for GLYCAM protein residues (NLN-to-neighbor)
- Loads `glycam-hydrogens.xml`, calls `Modeller.addHydrogens()` to restore stripped H
- 4-stage minimization: k=1000→100→10→0 kJ/mol/nm² restraints on protein heavy atoms; glycans move freely
- Configurable via `--relax-stages k1:iter1,k2:iter2,...`

**GLYCAM protein residues:** NLN (glycosylated ASN, ND2 bonds to sugar), OLS (glycosylated SER), OLT (glycosylated THR). These have peptide backbone but OpenMM doesn't recognize them as amino acids — bonds must be added explicitly.

**`--gromacs DIR`** exports GROMACS topology via ACPYPE pipeline:
1. OpenMM parametrizes with AMBER14 + GLYCAM_06j-1 → System (no constraints, so ParmEd gets all bond types)
2. ParmEd `load_topology()` → Structure → saves AMBER prmtop/inpcrd
3. ACPYPE converts to GROMACS `.top`/`.gro` with per-pair 1-4 parameters via `[ pairs_nb ]` directive (solves mixed AMBER fudgeLJ=0.5 / GLYCAM fudgeLJ=1.0 scaling)
4. Output: `{stem}.top`, `{stem}.gro`, `posre_{stem}.itp` in target directory

### dvbfixer puppet
Strips PDB to backbone-only polyglycine model. Removes all non-ATOM lines, keeps only backbone atoms (N, CA, C, O, OXT), renames all residues to GLY. `--keep CHAIN:NUM` preserves specific residues intact (all atoms, original name) — accepts single, range (`A:100-110`), list (`A:100,105`), or mixed, repeatable. No dependencies beyond stdlib.

### dvbfixer zbs
Full pipeline: renumber → model → prepare → minimize → protonate → minimize. Interim files are deleted by default (`--keep-interim` to preserve). Two minimize passes: first keeps existing H (default) for good heavy-atom positions, then protonate assigns AMBER names (HIE/GLH/CYX etc.) via PROPKA, then second minimize uses `--rebuild-h` to strip and re-add H with correct protonation variants (e.g. HE2 for GLH). Final output has AMBER protonation names; use `dvbfixer rename` for canonical PDB names. Internally a second protonate re-applies AMBER names after minimize (OpenMM's `PDBFile.writeFile` reverts them). Water removed by default. `.dat` flows from model → prepare (merged) → minimize. Supports `--mutate` for point mutations (passed to prepare step).

## Architecture Notes

- Each subcommand module has `parse_args(argv=None)` and `main(argv=None)` — the `argv` parameter allows the CLI dispatcher to pass subcommand arguments.
- `cli.py` dispatches `sys.argv[2:]` to the appropriate module's `main()`.
- Entry point defined in `pyproject.toml`: `dvbfixer = "dvbfixer.cli:main"`.

## GROMACS Topology Notes

- `top.py` parses RTP files directly — no dependency on `pdb2gmx` or GROMACS installation.
- FF files bundled in `FF/amber99sb-ildn-lipid21.ff/` and `FF/charmm36_ljpme-jul2022.ff/`. Used at build time for RTP parsing; output `.top` is self-contained.
- Output `.top` inlines all FF parameters: `[ defaults ]` from `forcefield.itp`, `[ atomtypes ]` from `ffnonbonded.itp`, bonded params from `ffbonded.itp`, `[ cmaptypes ]`/`[ nonbond_params ]` for CHARMM, water moleculetype (rigid settles), ion moleculetypes. Chain moleculetypes also inlined. Only `posre_*.itp` remain as separate `#include` (inside `#ifdef POSRES`).
- `_read_ff_content(path)` strips all preprocessor directives (#include, #define, #ifdef, etc.) for clean inlining.
- `_write_water_topology(f, path)` extracts only the rigid (settles) version from water .itp files that have `#ifndef FLEXIBLE` / `#else` blocks.
- RTP `[ bondedtypes ]` header defines function types for bonds/angles/dihedrals/impropers.
- Bond entries with `-C` mean previous residue's C atom, `+N` means next residue's N atom.
- AMBER has explicit NXXX/CXXX terminal entries in RTP; CHARMM uses TDB patch files.
- CHARMM carb.rtp has 363 sugar residues, all with full hydroxyl groups. No inter-residue bonds. Glycosidic linkages handled by removing HO + charge redistribution + atom type change.
- `[ intermolecular_interactions ]` must come after `[ molecules ]` in `.top` file. `gmx solvate` appends to end of file, so interchain SS must be appended after solvation.
- PDB atom name format: columns 13-16 = atom name (4 chars), column 17 = altLoc (space), columns 18-20 = residue name. 4-char atom names (HE21) start at column 13; shorter names start at column 14 with leading space.

## GLYCAM Integration Notes

- GLYCAM uses 3-character residue codes encoding linkage position + sugar identity + anomeric config (e.g. 4YB = beta-GlcNAc linked at O4, VMB = beta-Man linked at O3+O6).
- GLYCAM protein residues: NLN (glycosylated ASN), OLS (glycosylated SER), OLT (glycosylated THR).
- OpenMM's PDBFile does NOT infer intra-residue bonds for GLYCAM residues — must be added from FF templates.
- OpenMM's PDBFile does NOT add peptide bonds involving GLYCAM protein residues — must be added explicitly.
- NLN template expects only HD21 (not HD22) — ND2's other bond goes to sugar (external bond).
- GLYCAM fragments from GLYCAM-Web include terminal atoms (OXT, H2, H3) that must be stripped when fragment is mid-chain.
- `Modeller.loadHydrogenDefinitions('glycam-hydrogens.xml')` required before `addHydrogens()` for GLYCAM residues.
- CYS involved in disulfide bonds must be renamed to CYX before OpenMM parametrization (detected from CONECT SG-SG bonds).

## PDB Format Notes

- Residue identity = `(resSeq, iCode)` not just resSeq. Column 26 (0-based) is the insertion code.
- `set_resid()` must write both the 4-char resSeq (cols 22-25) AND clear the iCode (col 26).
- Antibody structures use Kabat/Chothia numbering: insertion codes at CDR loops (52A, 82A-C, 100A-J). These are NOT chain breaks.
- GROMACS PDB output often has blank chain IDs and continuous numbering across chains.
- `PDBFile.writeFile(..., keepIds=True)` is required to preserve chain IDs when writing through OpenMM.
- Glycan residues (BGL, BMA, NAG, etc.) have C and N atoms but NOT peptide bonds — C->N distance detection must be restricted to `STANDARD_RESIDUES`.
- CONECT serial remapping between PDB sources (acceptor vs graft) must use atom identity (chain, resseq, atomname) not serial numbers, since serials from different sources collide.
