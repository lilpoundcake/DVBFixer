# dvbfixer

A suite of Python CLI tools for preparing PDB (Protein Data Bank) structural biology files. Handles common issues with PDB files from MD simulations and the PDB database: missing chain IDs, antibody insertion codes, missing loops/residues, loop rebuilding with Modeller, multi-template homology modeling, energy minimization with selective restraints, protonation state assignment, GROMACS topology generation, GLYCAM glycoprotein transplanting, small molecule parametrization (GAFF2), and glycan conformational clustering from MD trajectories.

## Installation

Create the environment and install:

```bash
micromamba create -f environment.yml
micromamba activate dvbfixer
pip install -e .
```

**Modeller license:** The `model` command requires Modeller, which needs a free academic license key. Register at https://salilab.org/modeller/registration.html, then set the key in `<env>/lib/modeller-10.8/modlib/modeller/config.py`.

After installation, `dvbfixer` is available as a CLI command:

```bash
dvbfixer <command> [options]
```

Or without activating the environment:

```bash
micromamba run -n dvbfixer dvbfixer <command> [options]
```

## Commands

### dvbfixer split — Empirical Chain Splitting

Splits chains in PDB or GRO files that lack chain IDs (e.g. GROMACS MD output). Assigns unique chain IDs (A-Z, a-z, 0-9), inserts TER records, and renumbers residues per chain. GRO files are converted to PDB via MDAnalysis, preserving all residue names including protonation variants (GLUP, ASPP, etc.). Water, ions, and buffer particles (BUF/BUFF) are removed before chain detection to prevent false breaks, then optionally re-appended.

Chain breaks are detected by three criteria, applied in priority order:

1. **Residue number backward jump** — residue sequence number decreases (insertion codes like 82->82A are handled correctly and NOT treated as breaks)
2. **C->N peptide bond distance** — distance exceeds 2.5 A (any residue with backbone C/N atoms)
3. **Nearest-atom gap** — minimum distance between any atoms of consecutive residues exceeds 15 A (fallback for sugars, ligands, ions that lack peptide bonds)

#### Usage

```bash
# Basic usage — writes input_split.pdb
dvbfixer split input.pdb

# GRO file input (output is always PDB)
dvbfixer split simulation.gro -v

# Verbose output showing detected chains
dvbfixer split input.pdb -v

# Custom output and cutoffs
dvbfixer split input.pdb -o output.pdb -d 3.0 -g 20.0

# Disable distance-based detection (use only residue numbering)
dvbfixer split input.pdb --no-distance

# Keep original residue numbers
dvbfixer split input.pdb --no-renumber
```

#### Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o`, `--output` | `<input>_split.pdb` | Output file path |
| `-d`, `--distance-cutoff` | 2.5 | C->N peptide bond cutoff (angstroms) |
| `-g`, `--gap-cutoff` | 15.0 | Nearest-atom gap cutoff for non-protein residues (angstroms) |
| `--no-distance` | off | Disable all distance-based detection |
| `--no-renumber` | off | Keep original residue numbers |
| `--keep-water` | off | Keep water and ions in output (removed by default) |
| `-v`, `--verbose` | off | Print detected chain info |

---

### dvbfixer renumber — SEQRES-Based Residue Renumbering

Renumbers residues by aligning ATOM records to the SEQRES section via subsequence matching. Removes insertion codes (e.g. Kabat/Chothia antibody numbering 100A-J -> sequential 105-114) while preserving correct gap positions for missing residues.

Updates **all** PDB sections that reference residue numbers:
- ATOM, HETATM, TER
- HELIX, SHEET, SSBOND, LINK, CISPEP
- HET, DBREF, SEQADV
- CONECT (atom serial remapping)
- REMARK 465, 500, 610

#### Usage

```bash
# Basic usage — writes input_renum.pdb
dvbfixer renumber input.pdb

# Verbose output showing alignment details and gaps
dvbfixer renumber input.pdb -v

# Custom output
dvbfixer renumber input.pdb -o renumbered.pdb
```

#### Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o`, `--output` | `<input>_renum.pdb` | Output file path |
| `--keep-water` | off | Keep water molecules (HOH, WAT, TIP3, SOL) — removed by default |
| `--rename` | off | Rename non-canonical residues (AMBER/CHARMM) to standard names before processing |
| `-v`, `--verbose` | off | Print alignment details and gap positions |

#### How It Works

1. Parses SEQRES records to get the full sequence per chain
2. Extracts unique (resSeq, iCode, resname) tuples from ATOM records
3. Aligns ATOM residues to SEQRES as a subsequence — each ATOM residue name is matched to the next occurrence in SEQRES
4. Assigns new sequential numbering based on SEQRES position (position 1 = first SEQRES residue)
5. Non-SEQRES residues (waters, ligands) are numbered sequentially after the last SEQRES position
6. Chains without SEQRES entries are renumbered sequentially from 1
7. All PDB sections are updated with the new numbering

---

### dvbfixer model — Loop/Gap Rebuilding with Modeller

Rebuilds missing loops and gaps using Modeller's LoopModel. Identifies missing regions by aligning ATOM records to the SEQRES sequence (or a user-provided FASTA), then runs Modeller's loop modeling with MD refinement to fill them.

Non-protein chains (glycans, ligands) are included in the Modeller pipeline via `env.io.hetatm=True` with `'.'` (BLK residue) entries in the target sequence, so Modeller preserves them through loop modeling. Original chain IDs and residue numbering are restored automatically.

Writes a `.dat` file recording all atoms in rebuilt gap residues. This is merged by `dvbfixer prepare` with its own additions, so that `dvbfixer minimize` applies appropriate restraints to all rebuilt regions. Water molecules are removed by default (`--keep-water` to preserve).

#### Usage

```bash
# Basic usage — writes input_model.pdb
dvbfixer model input.pdb -v

# Higher quality (more sampling, slower)
dvbfixer model input.pdb --num-models 2 --num-loops 4 --md-level slow -v

# Use FASTA instead of SEQRES for complete sequence
dvbfixer model input.pdb --fasta sequence.fasta -v

# Keep Modeller working directory for debugging
dvbfixer model input.pdb --keep-workdir -v
```

#### Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o`, `--output` | `<input>_model.pdb` | Output file path |
| `--fasta` | none | FASTA file with complete sequence(s) (alternative to SEQRES) |
| `-n`, `--num-models` | 1 | Number of initial models to generate |
| `--num-loops` | 2 | Number of loop refinement models per initial model |
| `--md-level` | fast | MD refinement level: none, fast, slow, very_slow, slow_large |
| `--no-terminal` | off | Do not model missing N/C terminal residues (only rebuild internal gaps) |
| `--keep-water` | off | Keep water molecules (HOH, WAT, TIP3, SOL) — removed by default |
| `--keep-workdir` | off | Keep Modeller temp directory |
| `-v`, `--verbose` | off | Print Modeller progress |

#### How It Works

1. Parses SEQRES (or FASTA) for the complete target sequence
2. Builds a target PIR with protein chains from SEQRES and non-protein chains as `'.'` (BLK residues)
3. Reads the full PDB with `env.io.hetatm=True` so non-protein atoms are included
4. Creates a PIR alignment between the structure (with gaps as `-`) and the target sequence using Modeller's `align2d`
5. Runs `LoopModel` to generate initial model(s) filling the gaps
6. Refines loop regions with configurable MD level
7. Selects the best model by lowest `molpdf` score
8. Restores original chain IDs (Modeller A,B,C,... -> original letters)
9. Restores original residue numbering using the alignment (template positions get original numbers, gap-filled positions get interpolated numbers)

---

### dvbfixer prepare — Structure Fixing with PDBFixer

Adds missing residues, missing heavy atoms, and hydrogens using PDBFixer. User-specified protonation variants (HIE/HID/HIP/ASH/GLH/CYX from input PDB or `--mutate`) are passed as explicit OpenMM variants for correct hydrogen placement, then preserved in the output PDB via text-based post-processing. Standard HIS residues use pH-based auto-detection. Writes a `.dat` file recording which atoms were added (including variant overrides for downstream minimize). Automatically merges upstream `.dat` (from `dvbfixer model`) if present. Supports point mutations via `--mutate` including protonation variants (e.g. `--mutate A:83:HIP`).

#### Usage

```bash
# Basic usage — writes input_prepared.pdb and input_prepared.dat
dvbfixer prepare input.pdb -v

# Custom output
dvbfixer prepare input.pdb -o fixed.pdb --dat fixed.dat

# Keep crystallographic waters and heterogens
dvbfixer prepare input.pdb --keep-water --keep-heterogens

# Apply point mutations
dvbfixer prepare input.pdb --mutate A:39:ALA --mutate B:100:GLY -v
```

#### Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o`, `--output` | `<input>_prepared.pdb` | Output PDB file |
| `--dat` | `<output>.dat` | Restraint data file path |
| `--ph` | 7.0 | pH for hydrogen addition |
| `--keep-water` | off | Keep crystallographic waters |
| `--keep-heterogens` | off | Keep all heterogens (ligands, ions) |
| `--mutate` | none | Mutate a residue: CHAIN:RESNUM:NEW_AA (can be used multiple times) |
| `--rename` | off | Rename non-canonical residues (AMBER/CHARMM) to standard names before processing |
| `-v`, `--verbose` | off | Print detailed progress |

#### The .dat File

The `.dat` file is a JSON file that tracks which atoms were added/rebuilt. The model step writes a `.dat` recording gap-filled residue atoms. The prepare step merges that upstream `.dat` with its own additions (missing atoms, hydrogens). The minimize step uses the merged `.dat` to apply selective restraints.

Structure:

```json
{
  "description": "...",
  "total_added": 142,
  "residue_summary": {
    "A/GLY105": {"heavy": 4, "hydrogen": 3},
    "A/ALA106": {"heavy": 6, "hydrogen": 5}
  },
  "added_atoms": [
    {"chain": "A", "resid": "105", "icode": "", "resname": "GLY", "atom": "N", "element": "N"},
    ...
  ]
}
```

You can edit the `added_atoms` list to change which atoms receive weak/no restraints during minimization. For example, remove entries to make those atoms "original" (strong restraints), or add entries to make existing atoms "new" (weak/free).

---

### dvbfixer pull — Bond Pulling

Pulls atoms together to form bonds (disulfide bridges, glycosidic bonds) using OpenMM partial minimization. Atoms within a configurable radius of bond endpoints move freely; the rest are frozen via mass=0. Auto-removes conflicting hydrogens (CYS HG for disulfides, ASN HD22 for glycosidic bonds). Validates bond geometry before and after pulling.

#### Usage

```bash
# Form a disulfide bond between CYS residues
dvbfixer pull input.pdb --bond A:22:SG:A:96:SG -v

# Multiple bonds
dvbfixer pull input.pdb --bond A:22:SG:A:96:SG --bond A:300:ND2:A:1301:C1 -v

# Custom radius and output
dvbfixer pull input.pdb --bond A:22:SG:A:96:SG --radius 8.0 -o output.pdb
```

#### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--bond` | (required) | Bond specification: CHAIN1:RES1:ATOM1:CHAIN2:RES2:ATOM2 (repeatable) |
| `-o`, `--output` | `<input>_pull.pdb` | Output file path |
| `--radius` | 6.0 | Radius around bond endpoints for free atoms (angstroms) |
| `--target` | auto | Target bond distance (angstroms, auto-detected from element pair) |
| `--rename` | off | Rename non-canonical residues before processing |
| `-v`, `--verbose` | off | Print detailed progress |

---

### dvbfixer minimize — Energy Minimization with OpenMM

Energy-minimizes a PDB structure with OpenMM using selective restraints. Reads a `.dat` file (from `dvbfixer model` + `dvbfixer prepare`) to apply different restraint strengths to original vs newly added atoms. All HIS residues are automatically renamed to explicit variants (HIE/HID/HIP) before any OpenMM operation. By default, keeps existing hydrogens from input; if any residues are missing H (e.g. from mutation), uses PDBFixer to add them. With `--rebuild-h`, strips existing H, runs PDBFixer to fix missing heavy atoms (e.g. mutated residues) and terminal atoms (OXT for truncated chains), then re-adds correct H via OpenMM. Detects AMBER protonation names (HIE/GLH/CYX etc.) from the raw PDB and passes them as `variants` to `addHydrogens`, ensuring correct protonation hydrogens are added (e.g. HE2 for GLH).

#### Three-Tier Restraint System

| Atom category | Force constant | Purpose |
|---------------|----------------|---------|
| Original heavy atoms | 100 kcal/mol/A^2 | Keep resolved structure in place |
| New backbone (N, CA, C, O, CB) | 5 kcal/mol/A^2 | Maintain reasonable loop geometry |
| New sidechain + all hydrogens | 0 (free) | Full relaxation |

Minimization runs in two phases:
1. Full restraints (1000 iterations)
2. Restraints reduced 10x (1000 iterations)

#### Usage

```bash
# Minimize with .dat restraint info
dvbfixer minimize input_prepared.pdb --dat input_prepared.dat -v
# Outputs: input_prepared_minimized.pdb

# Vacuum minimization (no solvent box)
dvbfixer minimize input.pdb --no-solvent

# Without .dat — all atoms get strong restraints
dvbfixer minimize input.pdb
```

#### Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o`, `--output` | `<input>_minimized.pdb` | Output minimized PDB |
| `--dat` | `<input>.dat` | Restraint data file from `dvbfixer prepare` |
| `--ph` | 7.0 | pH for hydrogen addition if needed |
| `--ff` | amber19/protein.ff19SB.xml amber19/tip3p.xml | Force field XML files |
| `--padding` | 1.0 | Solvent padding in nm |
| `--restraint-k` | 100.0 | Strong restraint constant (kcal/mol/A^2) |
| `--weak-k` | 5.0 | Weak restraint constant for new backbone (kcal/mol/A^2) |
| `--max-iter` | 1000 | Max minimization iterations per phase |
| `--rebuild-h` | off | Strip and re-add hydrogens via OpenMM (default: keep existing) |
| `--no-solvent` | off | Minimize in vacuum |
| `--platform` | auto | OpenMM platform (CPU, CUDA, OpenCL, Reference) |
| `--rename` | off | Rename non-canonical residues (AMBER/CHARMM) to standard names before processing |
| `-v`, `--verbose` | off | Print detailed progress |

#### Recommended Workflow

```bash
# 1. Fix structure
dvbfixer prepare input.pdb -v
# -> input_prepared.pdb + input_prepared.dat

# 2. Inspect in PyMOL/VMD/ChimeraX, edit .dat if needed

# 3. Minimize
dvbfixer minimize input_prepared.pdb --dat input_prepared.dat -v
# -> input_prepared_minimized.pdb
```

---

### dvbfixer protonate — PROPKA3 Protonation Assignment

Runs PROPKA3 to predict per-residue pKa values, then renames titratable residues to their correct protonation state at the target pH and adds the corresponding hydrogen atoms. Uses AMBER force field naming conventions (HID/HIE/HIP, ASH, GLH, CYM, CYX, LYN). Existing hydrogens are stripped and re-added by OpenMM based on the renamed residue templates.

#### Protonation Logic

| Residue | Condition | Renamed to | Description |
|---------|-----------|------------|-------------|
| HIS | pKa > pH | HIP | Doubly protonated (charged) |
| HIS | pKa < pH | HIE (default) | Neutral, Ne2 protonated |
| HIS | pKa < pH | HID (option) | Neutral, Nd1 protonated |
| ASP | pKa > pH | ASH | Protonated (neutral) |
| GLU | pKa > pH | GLH | Protonated (neutral) |
| CYS | pKa >= 90 | CYX | Disulfide bridge |
| CYS | pKa < pH | CYM | Deprotonated thiolate |
| LYS | pKa < pH | LYN | Neutral |

#### Usage

```bash
# Basic usage — writes input_prot.pdb
dvbfixer protonate input.pdb

# Custom pH
dvbfixer protonate input.pdb --ph 6.5

# Full pKa summary table
dvbfixer protonate input.pdb --summary

# Show non-standard protonation changes
dvbfixer protonate input.pdb -v

# Use HID as default neutral histidine tautomer
dvbfixer protonate input.pdb --his-default HID
```

#### Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o`, `--output` | `<input>_prot.pdb` | Output file path |
| `--ph` | 7.0 | Target pH |
| `--his-default` | HIE | Default neutral HIS tautomer (HIE or HID) |
| `--cys-disulfide-pka` | 90.0 | pKa threshold for CYS -> CYX assignment |
| `--no-hydrogens` | off | Only rename residues, do not add/fix hydrogen atoms |
| `--ff` | amber19/protein.ff19SB.xml amber19/tip3p.xml | Force field XML files for hydrogen addition |
| `--keep-water` | off | Keep water molecules (HOH, WAT, TIP3, SOL) — removed by default |
| `--summary` | off | Print full pKa table |
| `-v`, `--verbose` | off | Print non-standard protonation changes |

---

### dvbfixer rename — Canonicalize Residue Names

Renames non-canonical residues to standard PDB three-letter codes. Handles AMBER protonation names (HIE/HID/HIP→HIS, ASH→ASP, GLH→GLU, CYX/CYM→CYS, LYN→LYS), CHARMM names (HSD/HSE/HSP→HIS), and selenomethionine (MSE→MET). Text-based — does not modify coordinates or atoms.

#### Usage

```bash
# Basic usage — writes input_renamed.pdb
dvbfixer rename input.pdb

# Verbose output
dvbfixer rename input.pdb -v
```

#### Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o`, `--output` | `<input>_renamed.pdb` | Output file path |
| `-v`, `--verbose` | off | Print each rename |

Also available as `--rename` flag on `renumber`, `prepare`, `minimize`, and `pull` tools to canonicalize input before processing.

---

### dvbfixer top — GROMACS Topology Generation

Generates GROMACS topology files directly from PDB or GRO files by parsing force field RTP/ARN/R2B/TDB files in Python — no `pdb2gmx` or GROMACS installation required. GRO files are auto-converted via MDAnalysis. Supports AMBER99SB-ILDN and CHARMM36 force fields (bundled). Output is a modular set of `.itp` files: `ffparams.itp` (all FF parameters), `{chain}.itp` (each chain moleculetype), `water.itp`, `ions.itp`, `posre_*.itp` (position restraints), `interchain_ss.itp` (inter-chain SS bonds + protein-glycan bonds), and a compact `topol.top` with only `#include` directives — no external FF directory needed. Handles proteins, carbohydrates (CHARMM glycan topology with glycosidic bond detection and protein-glycan bonds), glycolipids (ceramide + sugar tree as single moleculetype from CHARMM-GUI output), small CGenFF molecules (ACET, ACEH — auto-detected via distance splitting), lipids, nucleic acids, and all other CHARMM molecule types (~2400+ residues). Water (SOL/HOH/WAT), ions, and buffer particles (BUF) in the input are auto-detected and added to `[ molecules ]`. Automatically splits chains with overlapping residue numbers (e.g. duplicate glycan trees from `transplant`).

#### Usage

```bash
# Basic usage with AMBER (default) — writes topol.top, posre_*.itp, conf.pdb
dvbfixer top input.pdb

# GRO file input (auto-converted via MDAnalysis)
dvbfixer top system.gro --ff charmm

# Buffer system (ACET/ACEH auto-detected as small molecules)
dvbfixer top buffer.gro --ff charmm

# CHARMM36 force field
dvbfixer top input.pdb --ff charmm

# Explicit disulfide bonds and HIS protonation
dvbfixer top input.pdb --ss A:22:A:96 --his A:64:HID

# Protonate all ASP/GLU/HIS (CHARMM: ASPP/GLUP/HSP)
dvbfixer top input.pdb --ff charmm --protonate all

# Protonate specific residues
dvbfixer top input.pdb --ff charmm --protonate A:66,A:46:GLUP

# Glycolipid from CHARMM-GUI (ceramide + sugar tree as one moleculetype)
dvbfixer top glycolipid_charmm.pdb --ff charmm -o glycolipid_top/

# ACPYPE mode: full AMBER14+GLYCAM parametrization via OpenMM
# Handles mixed 1-4 scaling for glycoproteins via [ pairs_nb ]
dvbfixer top input.pdb --acpype

# Custom output
dvbfixer top input.pdb -o my_topology.top --pdb my_conf.pdb
```

#### Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o`, `--output` | `topol.top` | Output .top file |
| `--ff` | `amber` | Force field: `amber` or `charmm` |
| `--ff-dir` | (bundled) | Custom force field directory |
| `--water` | `tip3p` | Water model: tip3p, spc, spce, tip4p |
| `--ignh` | off | Ignore hydrogens in input PDB |
| `--ss` | auto | Disulfide bond: CHAIN1:NUM1:CHAIN2:NUM2 (repeatable) |
| `--his` | auto | HIS protonation: CHAIN:NUM:STATE (HIE/HID/HIP, repeatable) |
| `--protonate` | off | Protonate residues: `all` for every ASP/GLU/HIS, or `CHAIN:NUM[:STATE],...` for specific. H placed via OpenMM Modeller with CHARMM/AMBER FF |
| `--merge` | off | Merge all chains into single moleculetype |
| `--pdb` | `conf.pdb` | Output PDB with topology-matched atom names |
| `--acpype` | off | Use ACPYPE pipeline (AMBER14+GLYCAM -> ParmEd -> GROMACS) with per-pair 1-4 scaling |
| `-v`, `--verbose` | off | Print detailed progress |

#### RTP-based mode (default)

Parses force field RTP files to build topology from bond graph: resolves inter-residue bonds (`-C`, `+N`), enumerates angles/dihedrals/pairs algorithmically, copies impropers and CMAP from templates, applies terminal patches. Intra-chain disulfide bonds (SG-SG between CYS2 residues) are added explicitly to the chain topology with all derived angles, dihedrals, and 1-4 pairs; inter-chain SS bonds go to `interchain_ss.itp`. For CHARMM36, loads all molecule-type RTP files (aminoacids, carb, lipid, na, cgenff, ethers, metals, silicates, solvent — ~2400+ residue types). Non-protein chains are built without terminal patches. Glycan trees are detected from C1-O distances and built with proper glycosidic bond handling (HO removal, charge redistribution, atom type changes). Output is modular: `ffparams.itp` (all FF parameters including atomtypes, bonded params, water model params), chain `.itp` files, `water.itp`, `ions.itp`, `posre_*.itp` (position restraints, `#ifdef POSRES`), and `interchain_ss.itp` (if needed). `topol.top` contains only `#include` directives, `[ system ]`, and `[ molecules ]`. Ions and buffer particles (BUF) in the input PDB are auto-detected and added to the `[ molecules ]` section.

#### ACPYPE mode (`--acpype`)

Uses OpenMM to parametrize with AMBER14 + GLYCAM_06j-1, converts via ParmEd to AMBER prmtop/inpcrd, then ACPYPE generates GROMACS `topol.top`/`.gro` with `[ pairs_nb ]` directive for per-pair 1-4 parameters. This solves the mixed 1-4 scaling problem (AMBER fudgeLJ=0.5 vs GLYCAM fudgeLJ=1.0) that GROMACS cannot express globally. Output includes position restraints (`#ifdef POSRES` / `#include "posre_{stem}.itp"` / `#endif`) and water/ion moleculetypes ready for `gmx solvate`/`genion`. Best for glycoprotein systems. Ignores `--ff`/`--water`/`--merge` flags.

#### Glycolipid support (CHARMM36)

Automatically detects glycolipids — ceramide residues (CER1, CER160, CER180, etc.) covalently bonded to a sugar tree — and builds them as a single GROMACS moleculetype. Designed for CHARMM-GUI PDB output which uses 4-character residue names (CER1, BGLC, BGAL, ANE5, AFUC, etc.).

**How it works:**
- Detects sugar O atom within 2.0 A of ceramide C1S (CHARMM-GUI removes ceramide O1/HO1 at the bond site)
- Builds ceramide from `lipid.rtp` + sugar tree from `carb.rtp` into one moleculetype
- Handles linkage charge redistribution: ceramide O1+HO1 charge goes to C1S, sugar HO1 charge goes to O1
- Sugar O1 at ceramide junction gets type OC301 (linear ether); sugar-sugar junctions get OC3C61 (cyclic ether)
- Sialic acid (Neu5Ac/ANE5AC) handled correctly: links via C2, O2 removed at linkage
- Auto-detects BGAL with N-acetyl atoms (N, HN, CT) and remaps to BGALNA

**Supported ceramides:** CER160 (d18:1/16:0), CER180, CER181, CER2, CER200, CER220, CER240, CER241, CER3E.

**Example** (Fucosyl-GM1 ganglioside from CHARMM-GUI):
```bash
dvbfixer top FucGM1_Charmm.pdb --ff charmm -o gmx_top/ -v
# Output: Glycolipid_ _1.itp (251 atoms, 256 bonds, charge -1.0)
# Residues: CER160 + BGLC + BGAL + BGALNA + AFUC + ANE5AC
```

**Note:** Glycolipid support is CHARMM36 only. AMBER/GLYCAM+Lipid21 cannot be combined for glycolipids (incompatible atom type namespaces, no cross-FF parameters).

---

### dvbfixer transplant — Molecule Transplanting

Transplants molecules from a graft PDB into an acceptor PDB. Designed for the GLYCAM glycoprotein workflow: extract glycosylation site residues from your protein, submit to GLYCAM-Web, then transplant the GLYCAM output (renamed protein residues + glycan trees) back into the full structure. Also works with CHARMM-GUI output — use simple transplant mode (`--donor` + `--select`) to copy glycan chains or other molecules from CHARMM-GUI PDB into your structure.

#### Usage

```bash
# Graft workflow: replace donor residues in acceptor with GLYCAM output
dvbfixer transplant acceptor.pdb --donor donor.pdb --graft glycam_output.pdb

# With Kabsch superposition (if structures are not pre-aligned)
dvbfixer transplant acceptor.pdb --donor donor.pdb --graft glycam_output.pdb --superpose

# With OpenMM relaxation (AMBER+GLYCAM, 4-stage minimization)
dvbfixer transplant acceptor.pdb --donor donor.pdb --graft glycam_output.pdb --relax

# Export GROMACS topology via ACPYPE
dvbfixer transplant acceptor.pdb --donor donor.pdb --graft glycam_output.pdb --gromacs gmx_output/

# Simple transplant: copy selected residues from donor to acceptor
dvbfixer transplant acceptor.pdb --donor donor.pdb --select A:NAG

# CHARMM-GUI: transplant glycan chains from CHARMM-GUI output
dvbfixer transplant protein.pdb --donor charmm_gui.pdb --select G,H --superpose
```

#### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--donor` | (required) | Donor PDB: original residues from acceptor (identifies replacement sites) |
| `--graft` | none | Graft PDB: modified donor + added molecules (e.g. GLYCAM output) |
| `--select` | none | Selection for simple transplant: chain IDs, residue names, or ranges |
| `--align` | none | Chain mapping for superposition: DONOR:ACCEPTOR (e.g. H:H, repeatable) |
| `--superpose` | off | Enable Kabsch superposition (auto-detect chain mapping) |
| `--relax` | off | Run OpenMM minimization with AMBER+GLYCAM after transplant |
| `--relax-stages` | `1000:5000,...` | Relaxation stages: k1:iter1,k2:iter2,... (k in kJ/mol/nm2) |
| `--gromacs` | none | Export GROMACS topology via ACPYPE to specified directory |
| `-o`, `--output` | `<acceptor>_transplant.pdb` | Output PDB |
| `-v`, `--verbose` | off | Print detailed progress |

#### Graft Workflow

1. **`--donor`**: Original protein residues extracted from acceptor (e.g. ASN307, SER308). Used to identify which acceptor residues to replace and provides CA atoms for alignment.
2. **`--graft`**: GLYCAM output with renamed residues (NLN/OLS/OLT) + glycan trees (4YB/VMB/UYB etc.). Preserves GLYCAM atom and residue names.
3. Optional `--superpose`: Kabsch alignment of donor to acceptor, same transform applied to graft.
4. Donor residues removed from acceptor, graft protein residues inserted at correct positions, glycans appended.
5. Non-protein residues with resseq backward jumps are split into separate chains (prevents duplicate residue numbers when multiple glycan trees share a graft chain).
6. CONECT records remapped via atom identity (chain, resseq, atomname).

#### Relaxation (`--relax`)

4-stage energy minimization with AMBER14 + GLYCAM_06j-1:
- Protein heavy atoms restrained; glycans move freely
- Stages: k=1000 -> 100 -> 10 -> 0 kJ/mol/nm2
- Preprocessing: CYS->CYX for disulfides, bond addition for GLYCAM residues, hydrogen re-addition

#### GROMACS Export (`--gromacs DIR`)

Same ACPYPE pipeline as `dvbfixer top --acpype`: OpenMM parametrization -> ParmEd -> ACPYPE with `[ pairs_nb ]` for mixed 1-4 scaling. Outputs `topol.top`, `.gro`, and `posre_*.itp` to the specified directory. The `.top` includes position restraints and water/ion moleculetypes.

---

### dvbfixer glycam — Convert to GLYCAM Nomenclature

Converts PDB glycan structures from standard PDB sugar names (BGC, GAL, NAG, etc.) to GLYCAM force field 3-character codes encoding `[linkage][sugar][anomer]`. Detects glycosidic bonds from CONECT records (or distance-based fallback), determines linkage patterns, renames residues and atoms to GLYCAM convention (hydroxyl H: HO3→H3O; N-acetyl: C7→C2N, O7→O2N, C8→CME, N2→N2, HN→H2N; methyl: HT1→H1M etc.), and optionally adds ROH cap at the reducing end. Also detects protein-linked glycans and renames ASN→NLN, SER→OLS, THR→OLT. Writes ATOM records for protein residues and HETATM for sugars. Text-based — no OpenMM dependency. Handles input from PDB, CHARMM-GUI, or `dvbfixer prepare`.

#### Usage

```bash
# Basic usage — writes input_glycam.pdb
dvbfixer glycam glycan.pdb -v

# Without ROH cap at reducing end
dvbfixer glycam glycan.pdb --no-roh

# Custom output
dvbfixer glycam glycan.pdb -o glycam_output.pdb -v
```

#### Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o`, `--output` | `<input>_glycam.pdb` | Output file path |
| `--no-roh` | off | Do not add ROH cap at the reducing end |
| `-v`, `--verbose` | off | Print conversion details |

#### GLYCAM Naming Convention

Each sugar residue gets a 3-character name: `[linkage][sugar][anomer]`

**Linkage code** (1st character): `0`=terminal, `2`-`9`=single position, `V`=O3+O6, `W`=O3+O4, `U`=O4+O6, `Z`=O2+O3, `X`=O2+O6, `Y`=O2+O4 (multi-linkage for branching sugars).

**Sugar code** (2nd character): `G`=glucose, `L`=galactose, `M`=mannose, `Y`=GlcNAc, `V`=GalNAc, `f`=fucose (lowercase=L-sugar), `S`=Neu5Ac, `X`=xylose, `R`=ribose, etc.

**Anomer code** (3rd character): `A`=alpha, `B`=beta.

#### Example

```
BGC(res1, child at O4)     -> 4GB
GAL(res2, children at O3+O4) -> WLB
SIA(res6, terminal)        -> 0SA
NGA(res3, child at O3)     -> 3VB
GAL(res4, child at O2)     -> 2LB
FUC(res5, terminal)        -> 0fA
+ ROH cap at reducing end
```

---

### dvbfixer cluster — Glycan Conformational Clustering

Clusters glycan conformations from GROMACS MD trajectories using glycosidic torsion angle RMSD — the gold-standard method used by the Glycan Fragment Database (GFDB) and CHARMM-GUI Glycan Modeler. Auto-detects glycosidic linkages from the topology, extracts phi/psi/omega torsion angles across all frames, and clusters with the GROMOS algorithm.

Supports both CHARMM36 and GLYCAM force field naming, including sialic acid (Neu5Ac) which links via C2 instead of C1. Two clustering modes: `global` (all torsions at once) and `per-linkage` (each linkage independently, then combined into compound states — default, better for capturing per-linkage conformational variation). Representative structures are medoids (real frames closest to circular mean), automatically aligned on the root sugar or protein attachment point.

#### Usage

```bash
# Basic usage — auto-detects linkages, per-linkage clustering
dvbfixer cluster topology.tpr trajectory.xtc -v

# With interactive plots
dvbfixer cluster topology.tpr trajectory.xtc --plot -v

# Global clustering mode (GFDB-style, all torsions at once)
dvbfixer cluster topology.tpr trajectory.xtc --mode global

# Custom cutoff and stride
dvbfixer cluster topology.tpr trajectory.xtc --cutoff 20 --stride 10

# PDB input (no .tpr needed)
dvbfixer cluster structure.pdb trajectory.xtc -o my_clusters --plot

# Separate PDB per cluster, no alignment
dvbfixer cluster topology.tpr trajectory.xtc --separate-pdb --no-align

# Align on specific residue
dvbfixer cluster topology.tpr trajectory.xtc --align-resid 5
```

#### Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o`, `--output` | trajectory stem | Output file prefix |
| `--cutoff` | 30.0 | RMSD cutoff in degrees for GROMOS clustering |
| `--mode` | `per-linkage` | `global` (all torsions) or `per-linkage` (each linkage independently) |
| `--stride` | 1 | Read every Nth frame |
| `--begin` | 0 | First frame to read (0-based) |
| `--end` | last | Last frame (exclusive) |
| `--align-resid` | auto | Residue ID to align representatives on (auto: protein attachment or root sugar) |
| `--no-align` | off | Disable alignment of representative PDBs |
| `--separate-pdb` | off | Write each cluster as separate PDB (default: multi-MODEL PDB) |
| `--select` | all | MDAnalysis selection for output PDB atoms |
| `--plot` | off | Generate interactive HTML plots (requires plotly) |
| `-v`, `--verbose` | off | Verbose output |

#### Output Files

| File | Description |
|------|-------------|
| `{prefix}_torsions.csv` | Torsion angles per frame |
| `{prefix}_clusters.csv` | Cluster assignment per frame |
| `{prefix}_summary.txt` | Human-readable summary with per-cluster average torsions |
| `{prefix}_summary.json` | Machine-readable JSON summary |
| `{prefix}_representatives.pdb` | Multi-MODEL PDB with cluster representative structures (aligned) |
| `{prefix}_rama_{linkage}.html` | Ramachandran scatter + free energy surface (interactive, with `--plot`) |
| `{prefix}_timeseries.html` | Torsion angle time series colored by cluster (with `--plot`) |
| `{prefix}_populations.html` | Cluster population bar chart (with `--plot`) |

#### Torsion Angle Definitions

Crystallographic convention (IUPAC):

| Angle | Standard hexose | Sialic acid (Neu5Ac) |
|-------|----------------|---------------------|
| **phi** | O5–C1–Ox–C'x | O6–C2–Ox–C'x |
| **psi** | C1–Ox–C'x–C'(x-1) | C2–Ox–C'x–C'(x-1) |
| **omega** | Ox–C'6–C'5–O'5 (1→6 only) | same |

---

### dvbfixer homology — Multi-Template Homology Modeling

Multi-template homology modeling with Modeller. Takes a target FASTA (multi-chain) and one or more template PDB files. Auto-aligns target to templates via pairwise `align2d` per chain (or `--salign` for structure-based). Each target chain is modeled independently against its best template chain, then assembled into a multi-chain PDB. Point mutations are handled naturally by the differing target sequence. Antibody mode (`--antibody`): uses ANARCI for Kabat/IMGT numbering, CDR detection, and auto-mapping of Fv/constant domains to different templates.

#### Usage

```bash
# Basic multi-template
dvbfixer homology target.fasta --template fab.pdb --template fullsize.pdb -v

# With pipeline (prepare + minimize)
dvbfixer homology target.fasta --template fab.pdb --template fullsize.pdb --minimize -v

# Antibody mode
dvbfixer homology target.fasta --template fab.pdb --template igg.pdb --antibody -v
```

---

### dvbfixer parametrize — GAFF2 Small Molecule Parametrization

Parametrizes small molecules with GAFF2 force field and AM1-BCC or RESP charges for GROMACS MD. Wraps the AmberTools pipeline: antechamber → parmchk2 → tleap → ParmEd. Output: standalone `.itp` + `.gro` + `posre.itp`.

#### Usage

```bash
# AM1-BCC (default, fast)
dvbfixer parametrize molecule.pdb -n MOL -v

# Acetate with charge -1
dvbfixer parametrize acetate.pdb -n ACET --net-charge -1

# RESP (requires Gaussian log)
dvbfixer parametrize molecule.pdb -n MOL -c resp --gaussian-log molecule.log

# Generate Gaussian input for RESP
dvbfixer parametrize molecule.pdb -n MOL -c resp --gen-gaussian
```

---

### dvbfixer puppet — Backbone-Only Polyglycine Model

Strips a PDB to a minimal backbone scaffold: removes all non-ATOM lines, removes all sidechain and hydrogen atoms (keeps only N, CA, C, O, OXT), and renames every residue to GLY. Useful for creating "puppet" models for backbone-level alignment, modeling templates, or visualization.

#### Usage

```bash
# Basic usage — writes input_puppet.pdb
dvbfixer puppet input.pdb

# Keep specific residues intact (all atoms, original name)
dvbfixer puppet input.pdb --keep A:307

# Keep a range and a list
dvbfixer puppet input.pdb --keep H:286-296 --keep K:307,309

# Custom output
dvbfixer puppet input.pdb -o backbone.pdb
```

#### Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o`, `--output` | `<input>_puppet.pdb` | Output file path |
| `--keep` | none | Keep residue(s) intact: CHAIN:NUM, CHAIN:START-END, or CHAIN:NUM1,NUM2,START-END (repeatable) |

---

### dvbfixer zbs — Full Pipeline

Runs the complete preparation workflow in one command: **renumber → model → prepare → minimize → protonate → minimize**. Intermediate files are cleaned up by default — use `--keep-interim` to preserve them. Two minimize passes ensure correct protonation: the first keeps existing hydrogens (default) to get good heavy-atom positions, then protonate assigns AMBER protonation names (HIE/GLH/CYX etc.) based on PROPKA pKa predictions, then the second minimize uses `--rebuild-h` to strip and re-add hydrogens matching the correct protonation state (e.g. HE2 for GLH). The final output has AMBER protonation names — use `dvbfixer rename` if you need canonical PDB names. The `.dat` file flows from model (gap atoms) through prepare (merged with PDBFixer additions) to minimize (selective restraints). Water is removed by default. Each step can be skipped individually.

#### Usage

```bash
# Full pipeline
dvbfixer zbs input.pdb -v

# Skip terminal modeling and use vacuum minimization
dvbfixer zbs input.pdb --no-terminal --no-solvent -v

# Skip model and minimize steps
dvbfixer zbs input.pdb --skip-model --skip-minimize -v

# Custom output and pH
dvbfixer zbs input.pdb -o output.pdb --ph 6.5 -v

# With point mutations
dvbfixer zbs input.pdb --mutate A:39:ALA --mutate B:100:GLY -v
```

#### Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o`, `--output` | `<input>_zbs.pdb` | Final output PDB file |
| `--ph` | 7.0 | pH for protonation and hydrogen addition |
| `--ff` | amber19/protein.ff19SB.xml amber19/tip3p.xml | Force field XML files |
| `--skip-renumber` | off | Skip renumber step |
| `--skip-model` | off | Skip model step |
| `--no-terminal` | off | Do not model N/C terminal residues |
| `--num-loops` | 2 | Number of loop models |
| `--md-level` | fast | Modeller MD refinement level |
| `--fasta` | none | FASTA file for model step |
| `--skip-prepare` | off | Skip prepare step |
| `--keep-heterogens` | off | Keep heterogens during prepare |
| `--mutate` | none | Mutate a residue during prepare: CHAIN:RESNUM:NEW_AA (repeatable) |
| `--skip-minimize` | off | Skip minimize step |
| `--no-solvent` | off | Minimize in vacuum |
| `--rebuild-h` | off | Strip and re-add hydrogens via OpenMM during minimization (default: keep existing) |
| `--restraint-k` | 100.0 | Restraint force constant |
| `--max-iter` | 1000 | Max minimization iterations per phase |
| `--platform` | auto | OpenMM platform |
| `--skip-protonate` | off | Skip protonate step |
| `--no-hydrogens` | off | Only rename residues, skip hydrogen addition |
| `--keep-water` | off | Keep water molecules (removed by default) |
| `--keep-interim` | off | Keep all intermediate files (default: only final output) |
| `-v`, `--verbose` | off | Print detailed progress for all steps |

---

## Typical Pipeline

### Quick: Full pipeline with `zbs`

```bash
# Run everything: renumber -> model -> prepare -> minimize -> protonate -> minimize
dvbfixer zbs 6DDV.pdb -v
# -> 6DDV_zbs.pdb

# Skip terminal modeling for structures with missing N/C termini
dvbfixer zbs 6DDV.pdb --no-terminal -v

# Skip steps as needed
dvbfixer zbs input.pdb --skip-model --skip-minimize -v
```

### Manual: Step by step

For a GROMACS MD output PDB file with no chain IDs:

```bash
# 1. Split chains
dvbfixer split md_output.pdb -v

# 2. Renumber (requires SEQRES)
dvbfixer renumber structure_with_seqres.pdb -v

# 3. Rebuild missing loops (writes .dat for restraints)
dvbfixer model renumbered.pdb -v

# 4. Fix missing atoms/residues (merges model .dat)
dvbfixer prepare modeled.pdb -v

# 5. Minimize (pass 1 — standard protonation)
dvbfixer minimize prepared.pdb -v

# 6. Set protonation state names (rename only)
dvbfixer protonate minimized.pdb --no-hydrogens -v

# 7. Minimize (pass 2 — detects AMBER names, adds correct H)
dvbfixer minimize protonated.pdb --no-solvent -v

# 8. Re-apply AMBER names (OpenMM reverts them)
dvbfixer protonate minimized2.pdb --no-hydrogens -v
```

For a PDB database file with Kabat antibody numbering and missing loops:

```bash
# All at once (recommended)
dvbfixer zbs 1HZH.pdb -v

# Or step by step with custom options
dvbfixer renumber 1HZH.pdb -v
dvbfixer model 1HZH_renum.pdb --num-loops 4 -v
dvbfixer prepare 1HZH_renum_model.pdb -v
dvbfixer minimize 1HZH_renum_model_prepared.pdb -v
dvbfixer protonate 1HZH_renum_model_prepared_minimized.pdb --no-hydrogens -v
dvbfixer minimize 1HZH_renum_model_prepared_minimized_prot.pdb -v
dvbfixer protonate 1HZH_renum_model_prepared_minimized_prot_minimized.pdb --no-hydrogens -v
```

### Glycoprotein preparation with GLYCAM

#### Option A: Direct glycam → acpype pipeline (for structures with glycans already placed)

```bash
# 1. Convert PDB glycan names + atom names to GLYCAM convention
dvbfixer glycam glycoprotein.pdb -o glycam.pdb -v
# -> glycam.pdb (residues: UYB, 4YB, VMB, 0fA, NLN etc.; atoms: C2N, N2, H3O etc.)

# 2. Generate GROMACS topology (AMBER14 + GLYCAM, handles mixed 1-4 scaling)
dvbfixer top glycam.pdb --acpype -v
# -> topol.top, glycam.gro, posre_glycam.itp
```

#### Option B: GLYCAM-Web workflow (for adding glycans from scratch)

```bash
# 1. Prepare protein structure
dvbfixer zbs antibody.pdb -v
# -> antibody_zbs.pdb

# 2. Extract glycosylation sites into donor.pdb (manually or with PyMOL)
# 3. Submit donor.pdb to GLYCAM-Web -> get glycam_output.pdb

# 4. Transplant GLYCAM output back into full structure
dvbfixer transplant antibody_zbs.pdb --donor donor.pdb --graft glycam_output.pdb -v
# -> antibody_zbs_transplant.pdb

# 5. Relax with AMBER+GLYCAM force fields
dvbfixer transplant antibody_zbs.pdb --donor donor.pdb --graft glycam_output.pdb --relax -v
# -> antibody_zbs_transplant.pdb + antibody_zbs_transplant_relaxed.pdb

# 6. Generate GROMACS topology (handles mixed AMBER/GLYCAM 1-4 scaling)
dvbfixer top antibody_zbs_transplant.pdb --acpype
# -> topol.top, .gro, posre_*.itp
```

### GROMACS topology generation

```bash
# RTP-based (fast, modular .itp output — no FF dir needed)
dvbfixer top input.pdb --ff amber       # proteins
dvbfixer top input.pdb --ff charmm      # proteins + glycans + lipids + NA + more

# Glycolipid from CHARMM-GUI (auto-detected, CHARMM36 only)
dvbfixer top glycolipid_charmm.pdb --ff charmm -o gmx_top/

# ACPYPE-based (proteins + GLYCAM glycans, handles mixed 1-4 scaling)
dvbfixer top input.pdb --acpype
```

### Glycan conformational analysis

```bash
# Cluster glycan conformations from MD trajectory
dvbfixer cluster topol.tpr md.xtc --plot -v
# -> md_representatives.pdb, md_summary.json, interactive HTML plots
```

### Homology modeling (antibody engineering)

```bash
# Combine Fv from Fab template + constant domains from IgG template
dvbfixer homology target.fasta --template fab.pdb --template igg.pdb --antibody -v

# With full pipeline
dvbfixer homology target.fasta --template fab.pdb --template igg.pdb --minimize -v
```

### Small molecule parametrization

```bash
# Parametrize a buffer component for GROMACS
dvbfixer parametrize acetate.pdb -n ACET --net-charge -1
# -> ACET.itp, ACET.gro, posre_ACET.itp
```

---

## Known Issues

- **N-terminal ASH/GLH in ACPYPE mode**: AMBER14 has no N/C-terminal protonated ASP/GLU templates (NASH/NGLH — never parameterized via RESP in any AMBER version). When `--acpype` encounters ASH or GLH at chain termini, it strips the protonation hydrogen (HD2/HE2) and uses the standard deprotonated template (NASP/NGLU). A `UserWarning` is emitted. Internal (non-terminal) ASH/GLH residues are preserved correctly.

- **Chain ID mismatch in .dat workflow**: The `.dat` file stores chain IDs from PDBFixer. If the prepared PDB is saved through a tool that reassigns chain IDs (PyMOL, VMD), the `.dat` entries won't match the new chain letters. Workaround: ensure chain IDs remain consistent between prepare and minimize steps, or manually edit the `.dat` file.

- **Hydrogen handling in minimize**: By default, existing hydrogens are kept. Use `--rebuild-h` to strip and re-add via OpenMM (needed when protonation state changes). When AMBER protonation names (GLH, HIE, CYX, etc.) are detected in the input PDB, they are passed as `variants` to `addHydrogens` to ensure correct protonation hydrogens.

- **OpenMM normalizes AMBER names**: `PDBFile` reader converts GLH→GLU, HIE→HIS, CYX→CYS. The minimize tool reads raw PDB text first to capture original names. `PDBFile.writeFile` also writes standard names, so a final protonate text-based rename is needed to restore AMBER names.

- **Pull valence checking**: The `pull` tool validates bonds before and after pulling. Pre-pull: checks valence (bond count vs element max), warns about unusual element pairs. Post-pull: checks convergence (distance vs target), bond length range, and steric clashes within the pulling residues. All checks are warnings only — they do not prevent the operation.

- **Glycoprotein minimization in `minimize`**: Non-protein residues (glycans, ligands) are automatically stripped before minimization and restored with original coordinates afterward. For full glycoprotein minimization with GLYCAM force field, use `transplant --relax` instead.

- **Mixed 1-4 scaling (AMBER+GLYCAM)**: AMBER uses fudgeLJ=0.5/fudgeQQ=0.8333, GLYCAM uses 1.0/1.0. GROMACS only supports one global value. The `--acpype` flag on `top` and `--gromacs` on `transplant` solve this via ACPYPE's `[ pairs_nb ]` directive with per-pair LJ/Coulomb parameters.

- **AMBER14 has no terminal protonated ASP/GLU**: AMBER14 lacks NASH/NGLH/CASH/CGLH templates (no RESP charges were ever computed for terminal protonated ASP/GLU — a 15+ year gap). Affects both `dvbfixer top --acpype` and `dvbfixer top --ff amber --protonate`. When ASH/GLH is requested at a terminus, the protonation H is dropped, the residue is converted to standard ASP/GLU (using the existing NASP/CASP/NGLU/CGLU templates), and a `UserWarning` is emitted. HIS variants (HID/HIE/HIP) are unaffected — terminal templates exist (NHIE/CHIE etc.). CHARMM is unaffected — it uses TDB patches that combine cleanly with ASPP/GLUP.

- **Modeller terminal alignment**: `align2d` can misplace terminal gaps (e.g. matching last template residue to last target residue). This is auto-corrected by `_fix_terminal_alignment` which forces gaps to the actual N/C termini.

- **HIS tautomer selection**: PROPKA only predicts the overall pKa, not which nitrogen is protonated. The `--his-default` flag sets a global default (HIE or HID). For accurate per-residue tautomer assignment, use tools like MolProbity's Reduce or Schrodinger's ProtAssign.
