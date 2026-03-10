# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

dvbfixer is a Python package providing CLI tools for preparing PDB (Protein Data Bank) structural biology files. Installed as a single `dvbfixer` command with subcommands: `split`, `renumber`, `model`, `pull`, `prepare`, `minimize`, `protonate`, `rename`, `zbs`.

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
  zbs.py          — full pipeline (renumber->model->prepare->minimize->protonate->minimize)
  ffutils.py      — shared force field utilities (residue sets, OpenFF setup)
  pdbutils.py     — shared PDB utilities (CONECT remapping, serial maps)
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
Splits chains in PDB files lacking chain IDs (e.g. GROMACS output). Three detection criteria:
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
Energy-minimizes with OpenMM using selective restraints from the `.dat` file. Three tiers: original heavy=strong (100), new backbone=weak (5), new sidechain+H=free. Two phases: full restraints then 10x reduced. Non-protein residues (glycans, ligands) are stripped before parametrization and restored with original coordinates afterward. Before adding hydrogens, strips existing H and runs PDBFixer to auto-fix missing heavy atoms and terminal atoms (OXT) — handles mutated residues with wrong sidechains and truncated chains. `--keep-hydrogens` skips this and uses existing hydrogens. Detects AMBER protonation names (HIE/GLH/CYX etc.) from raw PDB text before OpenMM normalizes them, and passes them as `variants` to `addHydrogens` so correct protonation hydrogens are added (e.g. HE2 for GLH).

**Important:** OpenMM's `PDBFile` reader normalizes AMBER names (GLH→GLU, HIE→HIS, CYX→CYS). The raw PDB must be read first with `_read_amber_renames()` to capture original names before loading with PDBFile.

**Known issue:** `.dat` stores chain IDs from PDBFixer. External tools may reassign chain IDs, breaking `.dat` matching.

### dvbfixer protonate
Runs PROPKA3 for pKa prediction, renames residues to AMBER protonation names (HID/HIE/HIP, ASH, GLH, CYM, CYX, LYN) based on target pH. `--no-hydrogens` for text-based rename only (used in pipeline). Default mode strips H and re-adds via OpenMM with `variants` parameter.

### dvbfixer rename
Text-based rename of non-canonical residue names to standard PDB names. Converts AMBER protonation (HIE/HID/HIP→HIS, ASH→ASP, GLH→GLU, CYX/CYM→CYS, LYN→LYS), CHARMM (HSD/HSE/HSP→HIS), and MSE→MET. Also available as `--rename` flag on `renumber`, `prepare`, `minimize`, and `pull` tools. Uses `CANONICAL_MAP` dict in `rename.py`.

### dvbfixer zbs
Full pipeline: renumber → model → prepare → minimize → protonate(rename) → minimize → protonate(re-apply). Interim files are deleted by default (`--keep-interim` to preserve). Two minimize passes: first with standard protonation for good heavy-atom positions, then protonate renames to AMBER names, then second minimize detects AMBER names from raw PDB and passes `variants` to `addHydrogens` for correct protonation hydrogens (e.g. HE2 for GLH). Final protonate re-applies AMBER names via text-based rename since OpenMM's `PDBFile.writeFile` reverts them. Water removed by default. `.dat` flows from model → prepare (merged) → minimize. Supports `--mutate` for point mutations (passed to prepare step).

## Architecture Notes

- Each subcommand module has `parse_args(argv=None)` and `main(argv=None)` — the `argv` parameter allows the CLI dispatcher to pass subcommand arguments.
- `cli.py` dispatches `sys.argv[2:]` to the appropriate module's `main()`.
- Entry point defined in `pyproject.toml`: `dvbfixer = "dvbfixer.cli:main"`.

## PDB Format Notes

- Residue identity = `(resSeq, iCode)` not just resSeq. Column 26 (0-based) is the insertion code.
- `set_resid()` must write both the 4-char resSeq (cols 22-25) AND clear the iCode (col 26).
- Antibody structures use Kabat/Chothia numbering: insertion codes at CDR loops (52A, 82A-C, 100A-J). These are NOT chain breaks.
- GROMACS PDB output often has blank chain IDs and continuous numbering across chains.
- `PDBFile.writeFile(..., keepIds=True)` is required to preserve chain IDs when writing through OpenMM.
- Glycan residues (BGL, BMA, NAG, etc.) have C and N atoms but NOT peptide bonds — C->N distance detection must be restricted to `STANDARD_RESIDUES`.
