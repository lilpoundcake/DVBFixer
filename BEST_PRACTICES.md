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

**Switching the same source to CHARMM36**: after step 1 (or any later
step with GLYCAM names), invert with `--to-charmm` and use the RTP
path on `top`:

```bash
dvbfixer glycam glycam.pdb --to-charmm -o charmm.pdb -v
# -> NAG/NDG/BMA/MAN/GAL/FUL/SIA + ASN (standard PDB sugar codes),
#    with standard PDB atom names. Linkage info preserved via CONECT.
dvbfixer top charmm.pdb --ff charmm -o charmm_topol.top -v
# -> CHARMM36 modular .itp files
```

The reverse `--to-charmm` produces standard 3-char PDB sugar names
that both CHARMM-GUI and `dvbfixer top --ff charmm` recognize natively.
One source structure → two FFs (AMBER+GLYCAM via `--acpype` OR CHARMM36
via `--ff charmm`).

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

For loop modelling with a multi-chain FASTA (e.g. antibody Fv with heavy + light chains, or full IgG):

```bash
dvbfixer model input.pdb --fasta target.fasta -o modeled.pdb -v
```

The FASTA must use chain-ID headers — the mapping is by chain ID, not file order, so any order in the FASTA file works:

```
>chain_H
EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVA...
>chain_L
DIQMTQSPSSLSASVGDRVSITCRASQDVNTAVAWYQQKPGKAPKLLIY...
```

PDB-style headers like `>4HKZ_H` are also accepted (the trailing `_X` is treated as the chain ID).

**Antibody numbering** with `renumber --scheme`. Five schemes supported, all local — no web service. ANARCI handles V-domains; bundled human IgG1/Cκ/Cλ EU references handle C-domains via Needleman-Wunsch.

```bash
# Kabat numbering — CDR insertions at H100A/B/C etc.
dvbfixer renumber antibody.pdb --scheme kabat -v

# EU numbering across the full IgG (V 1-113 then 118-447)
dvbfixer renumber igg.pdb --scheme eu -v

# Mix per-chain (default chain assignment + Kabat for HL)
dvbfixer renumber bispecific.pdb \
    --chain-scheme H:kabat --chain-scheme L:kabat -v
```

Partial chains (Fc-only with CH2+CH3, Fv-only with just VH, hinge-truncated constructs) are handled automatically — the V-detection and C-alignment run independently. Non-antibody chains fall back to the SEQRES path.

When the V-scheme extends past EU position 117 (IMGT V ends at 128, Martin/Aho end at 113-149), the EU C-domain numbering is shifted forward by `(max_V_resseq + 5 - first_C_EU)` to avoid number collisions, and a warning is printed. If you need fully monotonic numbering without the shift, stick with `--scheme kabat`, `--scheme chothia`, or `--scheme eu`.

## Workflow 7 — Residue deletion mutations

Delete residues from a structure (e.g. for IgG truncation, hinge engineering, glycosite knockouts) via `prepare --mutate X:N:del`. The tool automatically handles three edge cases:

```bash
# Truncate the C-terminal lysine of an IgG heavy chain
dvbfixer prepare igg.pdb --mutate H:447:del -o trunc.pdb -v

# Knock out an N-glycosylation site — attached glycan tree disappears too
dvbfixer prepare glycoprotein.pdb --mutate A:297:del -o no_glycan.pdb -v

# Delete a CYS that forms a disulfide — partner is auto-reduced (CYX→CYS,
# HG regenerated)
dvbfixer prepare antibody.pdb --mutate H:22:del -o no_ss.pdb -v

# Multiple deletions (consecutive treated as one contiguous gap)
dvbfixer prepare igg.pdb --mutate H:446:del --mutate H:447:del -v

# Mix substitution and deletion
dvbfixer prepare in.pdb --mutate A:39:ALA --mutate H:446:del -v

# Insertion-code residues
dvbfixer prepare antibody.pdb --mutate H:100A:del -v
```

The `.dat` file records each deletion in a `removed_residues` field (chain, resid, icode, resname, gap_type, gap_distance_A, linked_glycan_residues, disulfide_partner_repaired). Downstream `minimize` / `model` reads this so it doesn't try to match against an upstream `.dat` that still listed the residue.

If a deletion produces an internal gap with `prev.C → next.N` distance > 5 Å, a warning is printed — consider running `dvbfixer pull` to close the gap or `dvbfixer model` to graft a replacement loop.

## Workflow 8 — Multi-state PDB / MD trajectory chain assignment

`dvbfixer split` now handles multi-MODEL PDBs (NMR ensembles, GROMACS trajectory exports with MODEL records) as one complex sampled at multiple states. The same chain IDs are reused in every MODEL (A B C in every MODEL, not A B C / D E F / G H I as a naive walk would produce):

```bash
dvbfixer split trajectory.pdb -v
# Detected 50 MODEL records, all with identical chain layout
# — reusing chain IDs across MODELs.
```

If the MODELs differ structurally (e.g. one MODEL has a ligand the others lack), the tool falls back to per-MODEL independent chain IDs and emits a warning.

## Workflow 9 — Glycan conformational clustering from MD

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
  ASN/SER/THR → NLN/OLS/OLT automatically **only when the bonded sugar is
  GLYCAM-named** (UYB/4YB/VMB/...). For PDB-named (NAG/NDG/BMA/...) or
  CHARMM-named (BGLC/BMAN/...) sugars, ASN keeps its standard name. NLN
  is a GLYCAM-specific convention; using it with non-GLYCAM sugars
  confuses downstream tools.

- **Glycosylation detection** is FF-agnostic in `prepare`: CONECT records
  + distance-based fallback (ASN ND2 / SER OG / THR OG1 within 2.0 Å of
  a sugar anomeric C). Works on inputs that have no CONECT for glycosidic
  bonds (some crystal PDBs, CHARMM-GUI output for certain residues). The
  extra HD22 on glycosylated ND2 is always removed, regardless of FF.

- **HETATM-tagged NLN/OLS/OLT in input**: `prepare` preprocesses such
  inputs by rewriting `HETATM` → `ATOM  ` for protein/GLYCAM residues
  before PDBFixer reads them. Otherwise OpenMM treats them as ligands
  and fails to infer peptide bonds to neighbours. Also drops spurious
  TER records between same-chain amino-acid residues (a TER forces a
  chain split). Both edits are no-ops on clean inputs.

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

- **Modeller BLK alignment error**: occurs when a chain ID appears in two
  disjoint file segments in the input PDB (e.g. chain A protein, then
  chain B sugar, then more chain A). Modeller segments by file-block and
  the PIR alignment fails. `dvbfixer model` now auto-reorders the input
  so every chain ID is contiguous before Modeller sees it — the error
  should be rare. If it still appears, check the input PDB for duplicated
  chain IDs in non-contiguous segments.

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
