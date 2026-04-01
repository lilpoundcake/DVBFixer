# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

dvbfixer is a Python package providing CLI tools for preparing PDB (Protein Data Bank) structural biology files. Installed as a single `dvbfixer` command with subcommands: `split`, `renumber`, `model`, `pull`, `prepare`, `minimize`, `protonate`, `rename`, `top`, `transplant`, `puppet`, `glycam`, `cluster`, `parametrize`, `homology`, `zbs`.

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
  glycam.py       — convert PDB glycan nomenclature to GLYCAM FF naming
  cluster.py      — glycan conformational clustering from MD trajectories
  parametrize.py  — GAFF2 small molecule parametrization (antechamber+tleap+ParmEd)
  homology.py     — multi-template homology modeling with Modeller (antibody mode with ANARCI)
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

Key packages: numpy, OpenMM 8.4, PDBFixer 1.12, scipy, PROPKA 3.5, Modeller 10.8, BioPython, MDAnalysis, plotly, ANARCI, HMMER.

**Modeller requires a license key** from https://salilab.org/modeller/registration.html (free for academics). Set in `<env>/lib/modeller-10.8/modlib/modeller/config.py`.

## Subcommands

### dvbfixer split
Splits chains in PDB or GRO files lacking chain IDs (e.g. GROMACS output). GRO files are converted to PDB via MDAnalysis (preserves all residue names including protonation variants like GLUP, ASPP). Water, ions, and buffer particles (BUF/BUFF) are stripped before chain detection to prevent false breaks (`--keep-water` re-appends them). Three detection criteria:
1. Residue number backward jump (insertion codes handled — equal resSeq with different iCode is NOT a break)
2. C->N peptide bond distance > 2.5 A (any residue with backbone C/N atoms — no name-based filtering)
3. Nearest-atom gap > 15 A (fallback for residues lacking C/N backbone atoms: sugars, ligands)

### dvbfixer renumber
Renumbers residues by aligning ATOM records to SEQRES via subsequence matching. Removes insertion codes (e.g. Kabat 100A-J -> sequential). Updates **all** PDB sections: ATOM, HETATM, TER, HELIX, SHEET, SSBOND, LINK, CISPEP, HET, DBREF, SEQADV, CONECT, REMARK 465/500/610. Each section has specific column positions — see the `update_*` functions and `remap_resid()` helper.

### dvbfixer model
Rebuilds missing loops/gaps using Modeller's LoopModel. Takes SEQRES (or --fasta) as complete sequence, aligns to ATOM records via `align2d`, runs loop modeling with configurable MD refinement. Non-protein chains (glycans, ligands) are included in the Modeller pipeline via `env.io.hetatm=True` with `'.'` (BLK) entries in the target PIR sequence so Modeller keeps them through loop modeling. Post-processing restores: (1) original chain IDs from Modeller's A,B,C,..., (2) original residue numbering with insertion codes using the alignment (template positions get original (resSeq, iCode), gaps get interpolated numbers). `--no-terminal` trims N/C terminal missing residues from the target sequence. Terminal alignment is auto-fixed after `align2d` to prevent misplaced terminal gaps. Writes a `.dat` file recording all atoms in rebuilt gap residues — prepare merges this with its own additions. Water removed by default (`--keep-water` to preserve). Protonation variant names (HIE/HID/HIP/ASH/GLH etc.) are renamed to standard (HIS/ASP/GLU) before Modeller reads the PDB (Modeller only knows standard names), then restored in the output. When no gaps detected, copies input to output without running Modeller.

### dvbfixer pull
Pulls atoms together to form bonds (disulfide bridges, glycosidic bonds) using OpenMM partial minimization. Supports multiple `--bond` specifications. Protein-protein bonds use `CustomBondForce`; protein-HETATM bonds use `CustomExternalForce` toward the fixed HETATM position. Atoms within `--radius` of bond endpoints are free to move (mass=0 freezing for the rest). Auto-removes conflicting hydrogens (CYS HG for disulfides, ASN HD22 for glycosidic bonds). Pre-pull validation checks valence (bond count vs `MAX_BONDS`) and bond type reasonableness for the pulling residues. Post-pull validation checks convergence (distance vs target), bond length range, and steric clashes within pulling residues.

### dvbfixer prepare
Runs PDBFixer to add missing residues, heavy atoms, and hydrogens. Strips H before PDBFixer processing to ensure clean template matching (H re-added by `addMissingHydrogens`). All HIS residues are renamed to explicit variants (HID/HIE/HIP) before hydrogen addition to bypass OpenMM's unreliable auto-detection — default is HIE, detection uses HD1/HE2 if present, falls back to heavy atom check (ND1/NE2). Writes a `.dat` file (JSON) recording which atoms were added. Merges upstream `.dat` (from model step) if present — atoms from gap rebuilding carry through to minimize. Detects glycosylated residues from CONECT records and removes extra hydrogens (e.g. ASN HD22 on glycosylated ND2). Supports `--mutate CHAIN:RESNUM:NEW_AA` for point mutations (can be used multiple times). `--rename` canonicalizes non-standard residue names before processing. User's protonation variant names (HIE/HID/HIP/ASH/GLH/CYX etc.) from input PDB or `--mutate` are preserved in output PDB and saved as `variant_overrides` in `.dat` for downstream tools.

### dvbfixer minimize
Energy-minimizes with OpenMM using selective restraints from the `.dat` file. Three tiers: original heavy=strong (100), new backbone=weak (5), new sidechain+H=free. Two phases: full restraints then 10x reduced. Non-protein residues (glycans, ligands) are stripped before parametrization and restored with original coordinates afterward. All HIS residues renamed to explicit variants (HIE/HID/HIP) before any OpenMM operation. Default: keeps existing hydrogens from input; if any residues are missing H (e.g. from mutation), uses PDBFixer to add them (OpenMM's `addHydrogens` can't handle residues with no H). `--rebuild-h` strips existing H, runs PDBFixer to fix missing heavy atoms and terminal atoms (OXT), then re-adds correct H via OpenMM. Detects AMBER protonation names (HIE/GLH/CYX etc.) from raw PDB text before OpenMM normalizes them, and passes them as `variants` to `addHydrogens` so correct protonation hydrogens are added (e.g. HE2 for GLH). Also reads `variant_overrides` from `.dat` file (saved by prepare) to recover protonation info even when output PDB has standard names.

**Important:** OpenMM's `PDBFile` reader normalizes AMBER names (GLH→GLU, HIE→HIS, CYX→CYS). The raw PDB must be read first with `_read_amber_renames()` to capture original names before loading with PDBFile.

**Important:** OpenMM's `PDBFile` reader normalizes AMBER names (GLH→GLU, HIE→HIS, CYX→CYS). The raw PDB must be read first with `_read_amber_renames()` to capture original names before loading with PDBFile.

**Known issue:** `.dat` stores chain IDs from PDBFixer. External tools may reassign chain IDs, breaking `.dat` matching.

### dvbfixer protonate
Runs PROPKA3 for pKa prediction, renames residues to AMBER protonation names (HID/HIE/HIP, ASH, GLH, CYM, CYX, LYN) based on target pH. `--no-hydrogens` for text-based rename only (used in pipeline). Default mode strips H and re-adds via OpenMM with `variants` parameter.

### dvbfixer rename
Text-based rename of non-canonical residue names to standard PDB names. Converts AMBER protonation (HIE/HID/HIP→HIS, ASH→ASP, GLH→GLU, CYX/CYM→CYS, LYN→LYS), CHARMM (HSD/HSE/HSP→HIS), and MSE→MET. Also available as `--rename` flag on `renumber`, `prepare`, `minimize`, and `pull` tools. Uses `CANONICAL_MAP` dict in `rename.py`.

### dvbfixer top
Generates GROMACS topology files directly from PDB or GRO files by parsing force field RTP/ARN/R2B/TDB files in Python (no `pdb2gmx`). GRO files are auto-converted to PDB via MDAnalysis. Supports AMBER99SB-ILDN and CHARMM36 force fields (bundled in `FF/` directory). Output is a modular set of `.itp` files with a compact `topol.top` containing only `#include` directives, `[ system ]`, and `[ molecules ]` — no external FF directory needed in the GROMACS working directory. Water (SOL/HOH/WAT/TIP3) and ions are auto-detected and counted for `[ molecules ]`. Small CGenFF molecules (ACET, ACEH, etc.) are auto-detected via distance-based chain splitting and built as independent moleculetypes.

**Topology generation algorithm:**
1. Resolve PDB residue names to RTP names via `CANONICAL_MAP` + R2B mapping
2. Map PDB atom names to FF atom names via ARN file
3. Build atom list from RTP `[ atoms ]` with types, charges, masses
4. Build bond graph from RTP `[ bonds ]`, resolving inter-residue refs (`-C`, `+N`)
5. Enumerate angles (all 3-atom paths), dihedrals (all 4-atom paths), 1-4 pairs
6. Copy impropers and CMAP (CHARMM) from RTP
7. Apply terminal patches (CHARMM: NH3+/COO- from TDB files; AMBER: separate NXXX/CXXX RTP entries)

**All CHARMM molecule types:** Loads all available RTP files at startup: `aminoacids.rtp`, `carb.rtp`, `lipid.rtp`, `na.rtp`, `cgenff.rtp`, `ethers.rtp`, `metals.rtp`, `silicates.rtp`, `solvent.rtp` (~2400+ residue types). Also loads corresponding `.r2b` files. Chain filter keeps any chain with residues in `STANDARD_AA`, `PDB_TO_GMX`, `builder.residues`, or `PDB_TO_CARB`. Non-protein chains are built without terminal patches or SS bond detection.

**Glycan support (CHARMM):** Parses `carb.rtp`, maps PDB sugar names (NAG, BMA, MAN, GAL, FUC, SIA, BGL, AFU, AMA, BGA) to CHARMM RTP names via `PDB_TO_CARB`. Also accepts CHARMM-GUI native 4-char names (BGLC, BGAL, AFUC, ANE5, etc.) via extended `PDB_TO_CARB` entries. Detects glycosidic bonds from C1-O distances < 2.0 A using parsed chain data (not raw PDB). Sialic acid (Neu5Ac) links via C2 anomeric carbon, not C1. Builds glycan trees via BFS. Handles charge redistribution (HO removal → charge to O) and atom type changes (OC311→OC3C61 at linkage sites). When linked O atoms are not in PDB (CHARMM-GUI removes them at glycosidic bond sites), their combined charge (O + HO) is redistributed to the anomeric carbon (C1 or C2). `_resolve_sugar_rtp()` auto-detects BGAL with N-acetyl atoms (N, HN, CT) and remaps to BGALNA. Inter-residue glycan angles use default ftype (5 for CHARMM, with Urey-Bradley) to match GROMACS parameter lookup. `_GLYCAN_LINKAGE_PARAMS` provides extra bond/angle/dihedral parameters for glycosidic linkage sites where OC311→OC3C61 creates atom type combos not in the standard CHARMM36 distribution (parameters by analogy with CC321D/CC321C variants). Also includes sialic acid C2 linkage params (CC3062-OC3C61 angles/dihedrals for Neu5Ac-galactose junctions).

**Glycolipid support (CHARMM):** Detects ceramide residues (`PDB_TO_LIPID` maps CHARMM-GUI names like CER1→CER160) bonded to sugar trees. Detection: sugar O atom within 2.0 A of ceramide C1S atom (CHARMM-GUI removes ceramide O1/HO1 at the glycosidic bond). `build_glycolipid_chain()` builds the entire glycolipid (ceramide + sugar tree) as a single moleculetype from `lipid.rtp` + `carb.rtp`. Linkage chemistry: ceramide O1+HO1 charge redistributed to C1S, sugar root O1 type changed to OC301 (linear ether, not OC3C61 which is sugar-sugar cyclic ether), sugar HO1 charge redistributed to O1. Inter-residue bond: ceramide C1S → sugar O1. Ceramide-sugar CTO2-OC301 bond/angle/dihedral parameters already exist in `ffbonded.itp` (no extras needed). `read_pdb_chains()` handles 4-char resnames (CER1, BGLC, etc.) that extend into PDB column 21. Chains containing ceramide+sugar are skipped in the main chain loop and built by the glycolipid builder. `_is_ceramide()` checks both PDB and RTP names. `CERAMIDE_RTP` set: CER160, CER180, CER181, CER2, CER200, CER220, CER240, CER241, CER3E.

**Chain splitting in `read_pdb_chains()`:** Detects resseq backward jumps within a chain (e.g. two glycan trees with same chain ID and overlapping residue numbers from `transplant`) and splits into separate sub-chains with generated chain IDs. Sugar-only chains are skipped in the main chain-building loop (handled by glycan detection).

**Protonation (`--protonate`):** `--protonate all` protonates ALL ASP→ASPP, GLU→GLUP, HIS→HSP (CHARMM) or ASP→ASH, GLU→GLH, HIS→HIP (AMBER). `--protonate CHAIN:NUM[:STATE],...` protonates specific residues (STATE defaults to protonated form based on FF). Missing protonation H atoms (e.g. HD2 for ASPP, HE2 for GLUP, HD1 for HSP) are added with proper geometry using OpenMM's `Modeller.addHydrogens(forcefield, variants=variants)` — same approach as the `protonate` tool. Uses CHARMM FF (`charmm36.xml`) or AMBER FF (`amber14-all.xml`) matching the `--ff` selection. OpenMM variant names (ASH, GLH, HIP, HID, HIE) work with both FFs. The full protein is loaded, H stripped, non-protein residues (glycans) removed, then `addHydrogens` places all H with correct geometry; only the specific protonation H atoms are extracted back into the chain data. Coordinates are near-optimal — should still be refined by energy minimization.

**Disulfide bonds:** Intra-chain SS bonds (SG-SG between CYS2 residues within the same chain) are added directly to the chain `.itp` `[ bonds ]` section by `build_chain()`. The SG-SG bond is not in the CYS2 RTP definition (it's an inter-residue special bond), so it's added explicitly from the detected SS pairs. The bond is added before the adjacency graph is built, so all derived angles (CB-SG-SG, SG-SG-CB), dihedrals, and 1-4 pairs are generated automatically. Inter-chain SS bonds are written to `interchain_ss.itp` with `[ intermolecular_interactions ]`. Protein-glycan bonds (ASN ND2 - NAG C1, r0=0.143 nm from CHARMM CC3162-NC2D1) also go in `interchain_ss.itp`. Auto-included in `topol.top` after `[ molecules ]` (GROMACS requires this directive after the molecules section). Warning: after `gmx solvate`/`genion`, the `#include` must be moved below SOL/ion entries to remain at the end.

**Output:** Modular set of files:
- `topol.top` — compact file with only `#include` directives, `[ system ]`, and `[ molecules ]`
- `ffparams.itp` — all FF parameters (`[ defaults ]`, `[ atomtypes ]`, `[ bondtypes ]`, `[ angletypes ]`, `[ dihedraltypes ]`, CMAP, NBFIX, glycan linkage params)
- `{chain_name}.itp` — each chain moleculetype in its own file (e.g. `Protein_chain_H.itp`, `Glycan_A_1001.itp`, `Glycolipid_ _1.itp`)
- `water.itp` — water moleculetype (rigid settles)
- `ions.itp` — ion moleculetypes
- `posre_*.itp` — position restraint files (used with `#ifdef POSRES`)
- `interchain_ss.itp` — inter-chain SS bonds + protein-glycan bonds with `[ intermolecular_interactions ]` (must stay at end of `topol.top`, after SOL/ions)
- `conf.pdb` — output PDB with FF atom names (includes ions/BUF particles from input)
FF directory is only used at build time for RTP parsing, not needed at runtime.

**Key writers:** `_write_moleculetype(f, chain_top, bonded_types)` writes a chain moleculetype section to a file. `write_top()` generates the modular output via `_read_ff_content()` (strips preprocessor directives), `_parse_defaults()` (extracts `[ defaults ]`), `_write_water_topology()` (extracts rigid settles version from water .itp), and `_dedup_atomtypes()` (removes duplicate water/ion atom type entries from #ifdef blocks, e.g. HT, OT with heavy/real mass variants).

**Key data structures:** `ChainTopology` dataclass holds atoms, bonds, pairs, angles, dihedrals, impropers, cmap per chain. `TopologyBuilder` class loads FF data, builds chains, writes output. `rtp_parser.py` provides `parse_rtp()`, `parse_r2b()`, `parse_arn()`, `parse_tdb()`, `parse_atomtypes()`.

**Atom name matching (`_match_atom_names`):** Multi-pass algorithm: (0) exact match, (0a) ARN reverse mapping (e.g. HN→H for CHARMM), (1) strip trailing digits, (2) common H renames, (3) numbered variants, (4) singleton numbered atoms (HG1→HG). Validates all RTP heavy atoms are present; warns on extra/missing atoms.

**`--acpype` mode:** Alternative to RTP-based topology. Uses `acpype_export.py` shared module: OpenMM (AMBER14+GLYCAM) → ParmEd → ACPYPE → GROMACS `topol.top`/`.gro` with `[ pairs_nb ]` for mixed 1-4 scaling. Output includes `#ifdef POSRES` / `#include "posre_{stem}.itp"` / `#endif` in the moleculetype section, and water/ion moleculetypes appended before `[ system ]`. Ignores `--ff`/`--water`/`--ignh`/`--merge` flags. Respects `--ss` for explicit disulfide bonds.

### dvbfixer transplant
Transplants molecules from a graft PDB into an acceptor PDB, with optional AMBER+GLYCAM energy minimization. Designed for the GLYCAM glycoprotein workflow: extract glycosylation site residues → submit to GLYCAM-Web → transplant back with glycans attached. Also works with CHARMM-GUI output via simple transplant mode (`--donor` + `--select`) to copy glycan chains or other molecules.

**Graft workflow** (`--donor` + `--graft`):
1. `--donor`: original protein residues extracted from acceptor (e.g. ASN307/SER308/THR309). Identifies which residues to replace and provides CA atoms for alignment.
2. `--graft`: GLYCAM output with renamed protein residues (NLN/OLS/OLT) + glycan trees (UYB/4YB/VMB etc.)
3. Kabsch superposition aligns donor→acceptor; same transform applied to graft
4. Donor residues removed from acceptor, graft protein residues inserted at correct positions, glycans appended
5. CONECT records remapped via (chain, resseq, atomname) identity matching (not serial numbers, which collide between sources)
6. `_renumber_graft()` detects resseq backward jumps in non-protein residues and assigns different chain IDs per segment (prevents duplicate residue numbers when multiple glycan trees share a graft chain)

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
4. Output: `topol.top`, `{stem}.gro`, `posre_{stem}.itp` in target directory. `topol.top` includes `#ifdef POSRES` / `#include "posre_{stem}.itp"` / `#endif` and water/ion moleculetypes.

### dvbfixer glycam
Converts PDB glycan nomenclature to GLYCAM force field naming. GLYCAM uses 3-character residue codes: `[linkage][sugar][anomer]`. Linkage detected from CONECT records (or C1-O distance < 2.0 Å fallback). Single linkage: `0`=terminal, `2`-`9`=position. Multi-linkage: `V`=O3+O6, `W`=O3+O4, `U`=O4+O6, etc. Sugar codes: `G`=Glc, `L`=Gal, `M`=Man, `Y`=GlcNAc, `V`=GalNAc, `f`=Fuc (lowercase=L-sugar), `S`=Neu5Ac. Anomer: `A`=alpha, `B`=beta. Handles sialic acid (anomeric C2 not C1), N-acetyl atom renames (C7→C2N, O7→O2N, C8→CME), ROH cap at reducing end, and protein-linked glycans (ASN→NLN, SER→OLS, THR→OLT). Text-based, no OpenMM dependency. H addition handled downstream by `transplant --relax`.

### dvbfixer homology
Multi-template homology modeling with Modeller. Takes a target FASTA (multi-chain) and one or more template PDB files. Auto-aligns target to templates via Modeller's `align2d` (or `salign` with `--salign`). Builds model with `automodel` or `LoopModel` using multiple `knowns`. Point mutations handled naturally by differing target sequence from templates. Post-processing restores chain IDs and residue numbering. Writes `.dat` file for downstream `prepare`/`minimize` restraints. `--prepare` and `--minimize` flags run the full pipeline automatically. Antibody mode (`--antibody`): uses ANARCI for Kabat/IMGT numbering, CDR detection, VH/VL/CH/CL domain classification, and auto-mapping of Fv from one template + constant domains from another. Dependencies: Modeller (required), ANARCI (for `--antibody` mode).

### dvbfixer parametrize
Parametrises small molecules with GAFF2 force field and AM1-BCC or RESP charges for GROMACS MD. Wraps the AmberTools pipeline: `antechamber` (atom types + charges) → `parmchk2` (missing parameter check) → `tleap` (AMBER topology) → ParmEd (AMBER→GROMACS conversion). Output: standalone `.itp` file (with `[ defaults ]`, `[ atomtypes ]`, `[ moleculetype ]` sections), `.gro` coordinates, and `posre_*.itp` position restraints. For RESP charges, user provides a Gaussian HF/6-31G* `.log` file (or uses `--gen-gaussian` to create the Gaussian input). AM1-BCC is the default (fast, no QM needed, ~95% of RESP accuracy). Supports PDB, MOL2, and SDF input formats.

### dvbfixer cluster
Clusters glycan conformations from MD trajectories using glycosidic torsion angle RMSD (GFDB method). Auto-detects glycosidic linkages from topology, extracts phi/psi/omega torsion angles (crystallographic convention: phi=O5-C1-Ox-C'x, psi=C1-Ox-C'x-C'(x-1), omega for 1→6 linkages), builds circular-RMSD distance matrix, runs GROMOS-style clustering. Handles both CHARMM36 and GLYCAM force field naming. Sialic acid uses C2 anomeric carbon and O6 ring oxygen. Two modes: `--mode global` (cluster all torsions simultaneously) and `--mode per-linkage` (cluster each linkage independently, combine into compound states — default, better at capturing per-linkage conformational variation). Representative structures are medoids (real frames closest to circular mean), aligned by Kabsch superposition on root sugar (auto-detected) or protein attachment point. Output: torsion CSV, cluster assignments CSV, JSON/text summary, representative PDBs (multi-MODEL or separate), interactive plotly HTML plots (Ramachandran + free energy surface, time series, population bar chart). Dependencies: MDAnalysis, numpy, plotly (for plots).

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
- FF files bundled in `FF/amber99sb-ildn-lipid21.ff/` and `FF/charmm36_ljpme-jul2022.ff/`. Used at build time for RTP parsing; output is modular `.itp` files — no external FF directory needed at runtime.
- Output is modular: `topol.top` has only `#include` directives + `[ system ]` + `[ molecules ]`. FF params in `ffparams.itp`, each chain in `{name}.itp`, water in `water.itp`, ions in `ions.itp`, position restraints in `posre_*.itp`, inter-chain SS in `interchain_ss.itp`.
- `_dedup_atomtypes()` removes duplicate atom type entries (e.g. HT, OT with heavy/real mass variants from #ifdef blocks).
- `_GLYCAN_LINKAGE_PARAMS` adds extra bond/angle/dihedral parameters for glycosidic linkage sites where OC311→OC3C61 creates atom type combos not in the standard CHARMM36 distribution (by analogy with CC321D/CC321C variants). Also includes sialic acid C2 linkage params (CC3062-OC3C61). Does NOT include ceramide-sugar CTO2-OC301 params (already in ffbonded.itp with proper multi-term dihedrals).
- `_read_ff_content(path)` strips all preprocessor directives (#include, #define, #ifdef, etc.) for clean inlining.
- `_write_water_topology(f, path)` extracts only the rigid (settles) version from water .itp files that have `#ifndef FLEXIBLE` / `#else` blocks.
- RTP `[ bondedtypes ]` header defines function types for bonds/angles/dihedrals/impropers.
- Bond entries with `-C` mean previous residue's C atom, `+N` means next residue's N atom.
- AMBER has explicit NXXX/CXXX terminal entries in RTP; CHARMM uses TDB patch files.
- CHARMM FF has ~2400+ residues across 9 RTP files: aminoacids (protein), carb (363 sugars), lipid (401), na (79 nucleic acids), cgenff (924 small molecules), ethers (25), metals (8), silicates (6), solvent/ions (77). All loaded at startup.
- Glycosidic linkages handled by removing HO + charge redistribution + atom type change. Linked O atoms not in PDB have their combined charge (O + HO) redistributed to the anomeric carbon (C1 or C2 for sialic acid).
- Glycolipids (ceramide + sugar tree) built as single moleculetype by `build_glycolipid_chain()`. Ceramide from `lipid.rtp`, sugars from `carb.rtp`. Ceramide-sugar bond: C1S—O1 (sugar O1 bridges). Ceramide O1+HO1 charge → C1S. Sugar O1 type: OC301 (not OC3C61). CTO2-OC301 params already in ffbonded.itp (multi-term dihedrals — must not be redefined). CHARMM-GUI 4-char resnames (CER1, BGLC, etc.) detected via `read_pdb_chains()` extended column parsing.
- `interchain_ss.itp` with `[ intermolecular_interactions ]` is auto-included in `topol.top` after `[ molecules ]` (GROMACS requires this directive after the molecules section). Contains both inter-chain SS bonds and protein-glycan bonds (ASN ND2 - NAG C1, r0=0.143 nm).
- Ions and buffer particles (BUF) are auto-detected in PDB by matching residue names against moleculetypes in `ions.itp`. Counted and added to `[ molecules ]` section. Not built as chain topologies — their moleculetypes are defined in `ions.itp`. BUF atomtype (dummy, no LJ) added to `ffnonbonded.itp`.
- Water molecules (SOL/HOH/WAT/TIP3) counted by atom count / 3 (not resseq dedup) to handle PDB resseq overflow in large systems (>10k residues wrap at 9999).
- Small CGenFF molecules (ACET, ACEH, ACEM, etc.) detected by distance-based chain splitting (`_split_chain_by_distance`, gap > 4 Å). Single-residue chains of known RTP types are counted and built as separate moleculetypes. Terminal patches are NOT applied to non-protein residues.
- GRO file support: auto-detected by `.gro` extension, converted to temp PDB via MDAnalysis (`_gro_to_pdb`). Original input path preserved for output directory and system name.
- 4-char resnames from GROMACS output (ACET, ACEH, TIP3) handled via `_KNOWN_4CHAR_RESNAMES` set in `read_pdb_chains()`.
- PDB serial numbers wrap at 100000 (`serial % 100000`) for large systems. Resseq in extra molecule lines renumbered sequentially.
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
- Glycan residues (BGL, BMA, NAG, etc.) have C and N atoms but NOT peptide bonds — `split_chains.py` detects breaks by backbone C/N atom presence (not residue name filtering).
- CONECT serial remapping between PDB sources (acceptor vs graft) must use atom identity (chain, resseq, atomname) not serial numbers, since serials from different sources collide.
