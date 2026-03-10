# dvbfixer

A suite of Python CLI tools for preparing PDB (Protein Data Bank) structural biology files. Handles common issues with PDB files from MD simulations and the PDB database: missing chain IDs, antibody insertion codes, missing loops/residues, loop rebuilding with Modeller, energy minimization with selective restraints, and protonation state assignment.

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

Splits chains in PDB files that lack chain IDs (e.g. GROMACS MD output). Assigns unique chain IDs (A-Z, a-z, 0-9), inserts TER records, and renumbers residues per chain.

Chain breaks are detected by three criteria, applied in priority order:

1. **Residue number backward jump** — residue sequence number decreases (insertion codes like 82->82A are handled correctly and NOT treated as breaks)
2. **C->N peptide bond distance** — distance exceeds 2.5 A (only for standard amino acid residues; glycan residues with C/N atoms are excluded)
3. **Nearest-atom gap** — minimum distance between any atoms of consecutive residues exceeds 15 A (fallback for sugars, ligands, ions that lack peptide bonds)

#### Usage

```bash
# Basic usage — writes input_split.pdb
dvbfixer split input.pdb

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
| `--keep-water` | off | Keep water molecules (HOH, WAT, TIP3, SOL) — removed by default |
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

Adds missing residues, missing heavy atoms, and hydrogens using PDBFixer. Writes a `.dat` file recording which atoms were added, for use by `dvbfixer minimize`. Automatically merges upstream `.dat` (from `dvbfixer model`) if present, so rebuilt gap atoms are tracked through the pipeline. Supports point mutations via `--mutate`.

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

### dvbfixer minimize — Energy Minimization with OpenMM

Energy-minimizes a PDB structure with OpenMM using selective restraints. Reads a `.dat` file (from `dvbfixer model` + `dvbfixer prepare`) to apply different restraint strengths to original vs newly added atoms. Before adding hydrogens, automatically strips existing H and runs PDBFixer to fix missing heavy atoms (e.g. mutated residues) and terminal atoms (OXT for truncated chains). Detects AMBER protonation names (HIE/GLH/CYX etc.) from the raw PDB and passes them as `variants` to `addHydrogens`, ensuring correct protonation hydrogens are added (e.g. HE2 for GLH).

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
| `--keep-hydrogens` | off | Use existing hydrogens from input (do not re-add) |
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

### dvbfixer zbs — Full Pipeline

Runs the complete preparation workflow in one command: **renumber → model → prepare → minimize → protonate → minimize → protonate**. Intermediate files are cleaned up by default — use `--keep-interim` to preserve them. Two minimize passes ensure correct protonation: the first minimizes with standard names to get good heavy-atom positions, protonate renames residues to AMBER names based on PROPKA pKa predictions, the second minimize detects AMBER names and re-adds hydrogens matching the correct protonation (e.g. HE2 for GLH), and the final protonate re-applies AMBER names (OpenMM reverts them on write). The `.dat` file flows from model (gap atoms) through prepare (merged with PDBFixer additions) to minimize (selective restraints). Water is removed by default. Each step can be skipped individually.

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
| `--keep-hydrogens` | off | Use existing hydrogens during minimization |
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

## Known Issues

- **Chain ID mismatch in .dat workflow**: The `.dat` file stores chain IDs from PDBFixer. If the prepared PDB is saved through a tool that reassigns chain IDs (PyMOL, VMD), the `.dat` entries won't match the new chain letters. Workaround: ensure chain IDs remain consistent between prepare and minimize steps, or manually edit the `.dat` file.

- **Hydrogen handling in minimize**: Hydrogens are stripped and re-added by OpenMM. Use `--keep-hydrogens` to preserve existing ones. When AMBER protonation names (GLH, HIE, CYX, etc.) are detected in the input PDB, they are passed as `variants` to `addHydrogens` to ensure correct protonation hydrogens.

- **OpenMM normalizes AMBER names**: `PDBFile` reader converts GLH→GLU, HIE→HIS, CYX→CYS. The minimize tool reads raw PDB text first to capture original names. `PDBFile.writeFile` also writes standard names, so a final protonate text-based rename is needed to restore AMBER names.

- **Pull valence checking**: The `pull` tool validates bonds before and after pulling. Pre-pull: checks valence (bond count vs element max), warns about unusual element pairs. Post-pull: checks convergence (distance vs target), bond length range, and steric clashes within the pulling residues. All checks are warnings only — they do not prevent the operation.

- **Glycoprotein minimization**: Non-protein residues (glycans, ligands) are automatically stripped before minimization and restored with original coordinates afterward. Full glycan minimization with force field is not currently supported.

- **Modeller terminal alignment**: `align2d` can misplace terminal gaps (e.g. matching last template residue to last target residue). This is auto-corrected by `_fix_terminal_alignment` which forces gaps to the actual N/C termini.

- **HIS tautomer selection**: PROPKA only predicts the overall pKa, not which nitrogen is protonated. The `--his-default` flag sets a global default (HIE or HID). For accurate per-residue tautomer assignment, use tools like MolProbity's Reduce or Schrodinger's ProtAssign.
