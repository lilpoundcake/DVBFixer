# Best practices for using dvbfixer

This guide describes the recommended end-to-end workflows for preparing
structures with dvbfixer. Each section is a self-contained recipe.

## Install

```bash
micromamba create -f environment.yml
micromamba activate dvbfixer
pip install -e .
```

After installation, `dvbfixer <command>` works. Modeller needs a free
academic license — set the key in `<env>/lib/modeller-10.8/modlib/modeller/config.py`.

## When to use which command

| Goal | Command |
|---|---|
| Fix missing chain IDs / split a GROMACS-style PDB | `dvbfixer split` |
| Fix antibody Kabat numbering / insertion codes | `dvbfixer renumber` |
| Rebuild missing loops/residues with Modeller | `dvbfixer model` |
| Add missing heavy atoms + protonate (PDBFixer) | `dvbfixer prepare` |
| Energy-minimize with selective restraints | `dvbfixer minimize` |
| Set protonation states via PROPKA3 pKa | `dvbfixer protonate` |
| Pull two atoms together (form SS, glycosidic) | `dvbfixer pull` |
| Convert glycan PDB → GLYCAM naming | `dvbfixer glycam` |
| Transplant residues between PDBs | `dvbfixer transplant` |
| Generate GROMACS topology | `dvbfixer top` |
| Parametrize a small molecule (GAFF2) | `dvbfixer parametrize` |
| Cluster glycan conformations from MD | `dvbfixer cluster` |
| Multi-template homology model | `dvbfixer homology` |
| Run the full prepare pipeline | `dvbfixer zbs` |

## Workflow 1 — Standard protein from PDB

The fastest path for a crystal structure with no glycans / unusual ligands:

```bash
dvbfixer zbs input.pdb -v
# renumber → model → prepare → minimize → protonate → minimize → protonate
# Final output: input_zbs.pdb (AMBER protonation names: HIE/GLH/CYX/...)
```

To get standard residue names back (HIS instead of HIE etc.):

```bash
dvbfixer rename input_zbs.pdb -o input_canonical.pdb
```

## Workflow 2 — Glycoprotein crystal structure (BioLuminate-style)

For quick visualisation or a one-shot energy-minimized PDB with PDB-named
glycans (NAG/FUC/GAL/MAN) covalently linked to ASN/SER/THR via CONECT
records. For production MD prefer **Workflow 3** (convert to GLYCAM
names first; gives proper AMBER+GLYCAM parametrization end-to-end).

```bash
# 1. Prepare adds H to protein AND glycans (RDKit-based, BioLuminate-style)
dvbfixer prepare crystal.pdb -o prep.pdb -v
# - CONECT records preserved in output (PyMOL/VMD show all bonds)
# - Glycosylated ASN renamed to NLN automatically
# - H placement respects external bonds (no extra H on C1 of N-linked NAG)

# 2. Energy-minimize protein with OpenMM, refine glycans with universal FF
dvbfixer minimize prep.pdb -o min.pdb --no-solvent \
    --obminimize-refine        # or --xtb-refine for higher quality
# - Protein minimized with AMBER14 (frozen heterogens during this pass)
# - Glycan heavy atoms + linkage refined by OpenBabel UFF/MMFF94 (auto-typing)
# - Protein anchor atoms frozen during refinement to preserve ASN-NAG bond
```

Result: every glycan H sits at ~1.1 Å from its parent C/O/N; glycosidic
bond geometry preserved (~1.45 Å for ASN ND2 – NAG C1).

`--obminimize-refine` is faster (seconds) and gives the cleanest bond
lengths. **Default FF is UFF** — MMFF94/MMFF94s mistypes the anomeric C
of N-linked sugars as sp2, giving incorrect 120° angles around the
glycosidic bond instead of the correct 109° sp3 tetrahedral.

`--xtb-refine` is slower (~minutes) and higher quality on paper, but
conda-forge currently ships xtb 6.7.1 which has a `$fix` bug that can
stretch the protein-glycan bond ~10 % — prefer obminimize until 6.8+
lands.

## Workflow 3 — GLYCAM-named glycoprotein (recommended for production MD)

When the input already uses GLYCAM 3-char sugar codes (UYB, 4YB, VMB,
0YA, 0fA, etc.) and glycoprotein residues (NLN/OLS/OLT) — typically the
output of `dvbfixer glycam` — full whole-system minimization works
without fallback. All four tools (`prepare`, `minimize`, `protonate`,
`top --acpype`) accept GLYCAM-named input end-to-end:

```bash
# 1. Convert PDB sugars to GLYCAM names (skip if already GLYCAM-named)
dvbfixer glycam crystal.pdb -o glycam.pdb -v
# - Detects glycosidic bonds from CONECT (or distance fallback)
# - Renames NAG→UYB/0YB, BMA→VMB, MAN→2MA, etc.
# - Renames ASN/SER/THR → NLN/OLS/OLT at glycosylation sites
# - Renames atoms (C7→C2N, O7→O2N, HO3→H3O, etc.)

# 2. Add hydrogens via AMBER14 + GLYCAM_06j-1 templates
dvbfixer prepare glycam.pdb -o prep.pdb -v
# - GLYCAM detection short-circuits RDKit/OpenBabel H polish (would
#   strip atoms with GLYCAM-specific names like C2N/O2N/CME)
# - add_glycam_bonds(positions=...) populates intra-residue + peptide
#   + sugar-sugar glycosidic bonds BEFORE addHydrogens
# - Filters NLN/OLS/OLT out of PDBFixer's nonstandard-residue warning

# 3. Energy-minimize the WHOLE system with AMBER14 + GLYCAM + SMIRNOFF
dvbfixer minimize prep.pdb -o min.pdb --no-solvent -v
# - ignoreExternalBonds=True for N-linked glycan junctions
# - Pre-solvent check uses residueTemplates to avoid CYM/CYX false
#   positives (no unnecessary fallback)
# - Heterogens free to relax; protein heavy atoms strongly restrained

# 4. Assign protonation states with PROPKA3
dvbfixer protonate min.pdb -o prot.pdb -v
# - NLN/OLS/OLT temporarily renamed to ASN/SER/THR for PROPKA; mapped
#   back afterwards
# - FF auto-switches to amber14+GLYCAM (ff19SB has no GLYCAM templates)
# - HID/HIE/HIP/CYX from input are preserved

# 5. Generate GROMACS topology (AMBER14+GLYCAM, mixed 1-4 scaling)
dvbfixer top prot.pdb --acpype -o gmx/ -v
# - Selective H stripping keeps HD21 (and other correctly-named H)
# - Defensive CYS template forcing handles inputs without SS CONECTs
# - Output: gmx/topol.top + gmx/prot.gro + gmx/posre_prot.itp
#   with [ pairs_nb ] for per-pair 1-4 scaling
```

GLYCAM names (NLN/UYB/4YB/VMB/0YB/...) and protonation variants
(HID/HIE/HIP/ASH/GLH/CYX/CYM/LYN) are preserved end-to-end.

## Workflow 4 — Custom protonation states

To set specific protonation states (e.g. from pKa calculations):

```bash
# Option A: use AMBER residue names in input PDB; prepare preserves them
# (input has HIE, ASH, GLH, CYX residues — they survive through the pipeline)

# Option B: --mutate flag with variant names
dvbfixer prepare input.pdb \
    --mutate A:83:HIP --mutate B:117:GLH --mutate A:34:ASH -v

# Option C: --protonate flag (in dvbfixer top for GROMACS topology)
dvbfixer top input.pdb --ff amber \
    --protonate H:6:GLH,L:1:ASH,H:101:HIP,H:208:HIE
```

For terminal residues, AMBER14 has no NASH/NGLH/CASH/CGLH templates —
dvbfixer detects this and emits a `UserWarning`, reverting to standard
ASP/GLU at the terminus.

## Workflow 5 — GROMACS topology generation

For MD-ready topology files:

```bash
# AMBER99SB-ILDN — fast RTP-based, no GROMACS install needed
dvbfixer top input.pdb --ff amber -o gmx/

# CHARMM36 — broader coverage (proteins, lipids, sugars, NA, ~2400 residues)
dvbfixer top input.pdb --ff charmm -o gmx/

# ACPYPE pipeline — AMBER14 + GLYCAM_06j-1 with mixed 1-4 scaling via [pairs_nb]
# Recommended for glycoproteins where you need GLYCAM sugar parameters.
dvbfixer top glycoprotein_glycam.pdb --acpype -o gmx/

# Output is modular .itp files + topol.top with only #include directives —
# no external FF directory needed at MD runtime.
```

For a glycoprotein from GLYCAM-Web through ACPYPE:

```bash
dvbfixer glycam crystal.pdb -o crystal_glycam.pdb       # rename to GLYCAM
dvbfixer top crystal_glycam.pdb --acpype -o gmx/
```

## Workflow 6 — Antibody-specific tasks

Multi-template homology model (e.g. Fab on one template + Fc on another):

```bash
dvbfixer homology target.fasta \
    --template fab.pdb --template igg1_fc.pdb \
    --antibody --minimize -v
# --antibody uses ANARCI for Kabat numbering + CDR detection
```

## Workflow 7 — Glycan conformational clustering from MD

After running MD on a glycoprotein:

```bash
dvbfixer cluster topol.tpr md.xtc --plot -v
# Auto-detects glycosidic linkages, extracts phi/psi/omega, runs GROMOS
# clustering per linkage, outputs Ramachandran HTML + PDB of representatives.
```

## Common gotchas

- **CONECT formatting in input PDB**: `prepare` canonicalizes malformed CONECT
  spacing automatically (e.g. `CONECT 5801 10293` where 5-digit serials don't
  fit the fixed-width columns). Without canonicalization OpenMM mis-parses
  these and produces spurious bonds.

- **NLN/OLS/OLT residues**: GLYCAM glycoprotein residues. `prepare` renames
  ASN/SER/THR → NLN/OLS/OLT automatically when a CONECT bond connects them
  to a sugar. Downstream tools (top, transplant, minimize) recognize these.

- **PDB sugar names with AMBER14**: NAG/FUC/etc. don't have AMBER14
  templates. `minimize` uses SMIRNOFF auto-parametrization for them via
  `create_forcefield_with_openff`. If that fails (e.g. unknown ligand without
  a SMILES entry in `KNOWN_GLYCAN_SMILES`), the tool auto-falls back to
  legacy strip-and-splice + `--obminimize-refine` for glycan geometry.

- **Disulfide bonds**: Detected from input CONECT records, with distance
  fallback (SG-SG within 2.5 Å) for inputs without CONECTs. Recognizes SG
  on CYS, CYX, and CYM. CYS residues in SS bonds are auto-renamed to CYX
  where needed. Inter-chain SS goes to `interchain_ss.itp` for the
  GROMACS topology output.

- **GLYCAM names survive end-to-end**: NLN/UYB/4YB/VMB/0YB/0fA/0LA/2MA
  pass through `prepare → minimize → protonate → top --acpype` without
  being renamed. `minimize.main()` snapshots NLN/OLS/OLT names from the
  raw input PDB at startup and restores them just before the final write
  — defensive belt-and-braces against the strip-and-splice fallback
  which internally renames them.

- **ATOM vs HETATM records**: OpenMM's `PDBFile.writeFile` defaults
  non-standard residue names (HID/HIE/HIP/ASH/GLH/CYX/CYM/LYN +
  NLN/OLS/OLT) to HETATM. All three of prepare/minimize/protonate
  post-process the output with `fix_atom_hetatm_records` to rewrite
  these back to ATOM. Don't expect HETATM lines for protein residues
  in the final output.

- **HD21 geometry preserved through `top --acpype`**: The `top --acpype`
  pipeline uses a per-residue allowlist `_GLYCAM_KEEP_H` for H atoms on
  NLN/OLS/OLT. Only disallowed atoms (HD22 on NLN, HG on OLS, HG1 on
  OLT, CHARMM HN/HT*) are stripped; HD21 and other correctly-named H
  survive into addHydrogens which preserves them. Result: the canonical
  trans-amide HD21 position from minimize (CG-ND2-HD21 ≈ 120°) is
  preserved in the output `.gro`.

- **Sugar-sugar bonds**: `add_glycam_bonds(positions=...)` detects
  glycosidic bonds by distance (anomeric C1, or C2 for sialic, within
  2.0 Å of a linkage O2/O3/O4/O6 on another sugar). Without this,
  GLYCAM templates for linkage-position sugars (e.g. 6LB declares
  O6 externally bonded) fail to match.

- **Terminal ASH/GLH (AMBER)**: AMBER14 has no NASH/CASH/NGLH/CGLH
  templates (15-year-old known gap — no RESP charges were ever computed).
  At terminals these are auto-reverted to standard ASP/GLU with the H
  stripped, plus a `UserWarning`.

- **Heterogen H drift during xtb/obminimize**: Fixed in current version —
  every new H atom RDKit places is bonded to its parent in the topology,
  so the CONECT record carries the bond into OpenBabel/xtb and they apply
  the proper bond-length restraint.

## Recommended pipeline for production MD

For a glycoprotein crystal structure → GROMACS-ready topology:

```bash
# 1. Convert to GLYCAM names (skip if already GLYCAM-named)
dvbfixer glycam   crystal.pdb -o glycam.pdb -v

# 2. Structure prep with AMBER14+GLYCAM templates
dvbfixer prepare  glycam.pdb -o prep.pdb -v

# 3. Energy minimize (full system, AMBER14+GLYCAM+SMIRNOFF)
dvbfixer minimize prep.pdb -o min.pdb --no-solvent -v

# 4. Assign protonation states via PROPKA
dvbfixer protonate min.pdb -o prot.pdb -v

# 5. GROMACS topology via ACPYPE (handles mixed AMBER+GLYCAM 1-4 scaling)
dvbfixer top      prot.pdb --acpype -o gmx_top/ -v

# 6. Standard GROMACS workflow with the generated topology
cd gmx_top/
gmx editconf -f prot.gro -o boxed.gro -bt cubic -d 1.0
gmx solvate -cp boxed.gro -p topol.top -o solv.gro
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -p topol.top -o ions.gro -neutral -conc 0.15
gmx grompp -f em.mdp -c ions.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em
# ... NVT, NPT, production
```

For a protein-only structure (no glycans): skip steps 1, drop `--acpype`
in step 5 (use `--ff amber` for the RTP path which produces modular
`.itp` files without needing the AMBER+GLYCAM topology merging).

Note for the ACPYPE pipeline: `topol.top` includes an `interchain_ss.itp`
`#include` at the end (with `[intermolecular_interactions]`). After
`gmx solvate`/`gmx genion`, move that include below the SOL/ion entries
in `topol.top` — GROMACS requires the directive to be last.
