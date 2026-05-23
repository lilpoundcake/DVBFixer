# dvbfixer architecture

## Module structure

```
src/dvbfixer/
├── cli.py              90 lines    — single entry point, dispatches to subcommand main()
├── __init__.py          3 lines    — package marker
│
│   STRUCTURE PREP PIPELINE (composable subcommands)
├── split_chains.py    356 lines    — empirical chain splitting (gap / dist / numbering)
├── renumber.py        446 lines    — SEQRES-based renumbering, removes insertion codes
├── model.py          1310 lines    — Modeller LoopModel loop/gap rebuilding
├── prepare.py        1344 lines    — PDBFixer wrapper + BioLuminate-style H placement
├── minimize.py       1200 lines    — OpenMM minimization + optional xtb/obminimize refine
├── protonate.py       387 lines    — PROPKA3 pKa-based protonation
├── pull.py            636 lines    — bond pulling via OpenMM mass=0 partial min
├── rename.py          102 lines    — text-based variant → canonical name
├── puppet.py           98 lines    — strip to backbone polyglycine
├── zbs.py             262 lines    — full pipeline (renumber→model→prepare→minimize→...)
│
│   FF / TOPOLOGY GENERATION
├── top.py            3507 lines    — RTP-based GROMACS topology (AMBER + CHARMM)
├── rtp_parser.py      261 lines    — parses GROMACS RTP/ARN/R2B/TDB/ATP files
├── acpype_export.py   760 lines    — ACPYPE-based GMX topology (OpenMM→ParmEd→ACPYPE)
├── ffutils.py         137 lines    — shared FF utilities + SMIRNOFF auto-parametrization
│
│   GLYCAN / SMALL-MOLECULE TOOLS
├── glycam.py          820 lines    — bidirectional PDB/CHARMM ↔ GLYCAM nomenclature converter
├── transplant.py      841 lines    — graft residues between PDBs (Kabsch align)
├── parametrize.py     406 lines    — GAFF2 small molecule (antechamber→tleap→ParmEd)
├── cluster.py        1177 lines    — glycosidic torsion clustering from MD trajectory
│
│   ANTIBODY / HOMOLOGY
├── homology.py        766 lines    — multi-template homology with ANARCI antibody mode
│
│   SHARED HELPERS
└── pdbutils.py         83 lines    — CONECT record remapping, atom serial maps
```

## Data flow

```
input.pdb ──split──→ split.pdb
            │            │
            │            └── renumber ──→ renumbered.pdb  (+ SEQRES → ATOM alignment)
            │                                  │
            │                                  └── model ──→ modeled.pdb + .dat (gap atoms)
            │                                                  │
            │                                                  └── prepare ──→ prepared.pdb + .dat (merged)
            │                                                                          │
            │                                                                          └── minimize ──→ minimized.pdb
            │
            └── For glycoproteins: glycam ─→ ... ─→ transplant ─→ minimize ─→ top
```

The `.dat` file is the structured handoff between pipeline stages. It is
JSON with a list of "added atoms" (chain, resid, icode, atom, element)
plus `variant_overrides` (protonation state). Each tool merges its
additions with upstream `.dat` so `minimize` can apply tiered restraints:

| Atom class | Source | Restraint |
|---|---|---|
| Original heavy | input crystal | strong (100 kcal/mol/Å²) |
| New backbone | model + prepare | weak (5 kcal/mol/Å²) |
| New sidechain + all H | prepare | free |
| Heterogens (recent change) | input or prepare | free (refined separately) |

## Three force-field paths

dvbfixer supports three independent paths to a GROMACS topology, chosen
by `dvbfixer top` flags:

### 1. RTP path (`--ff amber` or `--ff charmm`) — `top.py` + `rtp_parser.py`

Parses GROMACS `.rtp`/`.arn`/`.r2b`/`.tdb`/`.atp` files directly in Python
(no `pdb2gmx`, no GROMACS install). Builds bond graphs, enumerates angles
and dihedrals, applies terminal patches. CHARMM path loads ~2400 residue
templates across 9 RTP files (aminoacids, carb, lipid, na, cgenff,
ethers, metals, silicates, solvent).

Output: modular `.itp` files (`ffparams.itp`, per-chain `{chain}.itp`,
`water.itp`, `ions.itp`, `posre_*.itp`, `interchain_ss.itp`) with a
compact `topol.top` containing only `#include` directives. No FF dir
needed at runtime.

Strengths: fast, deterministic, no per-ligand setup. Strong glycan
support via CHARMM `carb.rtp`. Glycolipid support via auto-detection
of ceramide+sugar trees.

Weaknesses: only handles residues with RTP templates. Unknown organic
molecules fail.

### 2. ACPYPE path (`--acpype`) — `acpype_export.py`

OpenMM (AMBER14 + GLYCAM_06j-1) → ParmEd → ACPYPE → GROMACS. Solves
the mixed 1-4 scaling problem (AMBER fudgeLJ=0.5 vs GLYCAM fudgeLJ=1.0)
via ACPYPE's `[pairs_nb]` directive with per-pair LJ/Coulomb parameters.

Best for glycoprotein systems. Auto-detects SS bonds from CONECT,
reorders chains (protein first, glycan after), handles AMBER variant
names across PDBFile normalization, adds peptide bonds for GLYCAM
protein residues, uses `ignoreExternalBonds=True` + `residueTemplates`
for CYX disambiguation.

### 3. OpenFF auto-parametrization (used inside `minimize.py`) — `ffutils.create_forcefield_with_openff`

For OpenMM `createSystem` only (not for GROMACS topology output). Loads
AMBER14 + GLYCAM_06j-1 + registers `SMIRNOFFTemplateGenerator` for any
residue not in `PROTEIN_RESIDUES | SOLVENT_IONS` whose name appears in
`KNOWN_GLYCAN_SMILES`. Suppresses GLYCAM sugar templates when PDB-named
sugars are present (otherwise NAG fuzzy-matches to GLYCAM "QVA" with
mismatched atom set).

## Key abstractions

### `model.build_resnum_mapping` — three-stage resnum reconstruction

After Modeller produces `output.B9999*.pdb` with renumbered residues
(1..N) and reassigned chain IDs (A,B,C,...), `build_resnum_mapping`
reconstructs the original (resSeq, iCode) numbering per chain through a
three-stage pipeline:

1. `_find_seqres_offset_by_resseq` (K-finder) — deterministic resseq
   offset by letter-match tolerance. Primary path; fast and unambiguous
   when input resseq jumps line up with SEQRES missing positions.
   Tolerates up to 10% letter mismatches (mutations).
2. `_align_atoms_to_seqres` (NW fallback) — semi-global Needleman-Wunsch
   with affine gap penalties (open=-10, extend=-1, X-neutral, free end
   gaps on SEQRES side). Mutation-tolerant; handles chains where the
   deterministic offset doesn't exist.
3. align2d mask consumption (legacy) — original PIR-mask based path,
   kept as a final fallback for chains that defeat both deterministic
   paths (icode-bearing antibody Kabat numbering, etc.).

`_interpolate_gaps` is the reusable helper for all three stages: it
fills N-terminal, internal, C-terminal, and all-gap regions in one
place, given the flanking original (resSeq, iCode) tuples.

HETATM resseqs are preserved verbatim on the deterministic paths (no
longer renumbered sequentially after the protein), so ligand and glycan
numbering from the input PDB survives Modeller round-tripping.

### `model._reorder_chains_for_modeller` — pre-Modeller chain grouping

Modeller segments chains by file-block when reading a PDB. If a single
chain ID appears in two disjoint blocks (e.g. chain A protein + chain
B sugar + chain A more protein), Modeller treats the two A-blocks as
separate chains and the PIR alignment fails with a BLK alignment error
(or silently drops the second block). This preprocessing pass groups
every chain ID's ATOM/HETATM records into a single contiguous block in
the temp PDB Modeller reads. No-op on already-contiguous inputs.

### `model.parse_fasta` — chain-ID-keyed FASTA parsing

Returns `dict[chain_id, sequence]` (previously a list of tuples in
file order). Chain IDs are parsed from the FASTA header via priority
patterns: `chain_X`, `PDBID_X` (PDB-style `>4ABC_H`), then a final
bare `X` fallback. Chain order in the FASTA file is no longer
significant — the mapping is by chain ID, so users can list chains in
any order.

### `prepare.add_heterogen_h_via_rdkit` — BioLuminate-style H placement

Pipes the prepared PDB through RDKit:
1. `MolFromPDBBlock(proximityBonding=True)` — full molecular graph from
   CONECT + distance perception
2. Hardcoded C=O bond-order fixups for known PDB sugars (NAG/SIA)
3. `AddHs(addCoords=True)` — places H with proper 3D geometry
4. Post-filter: skip H on carbonyl atoms (per `_NO_H_ATOM_NAMES`)
5. Map RDKit atom indices → OpenMM atoms via positional lookup
6. Insert H atoms + **bond them to their parent atom**
7. Carry over RDKit-perceived intra-sugar bonds → OpenMM topology

The H-parent bond is critical: without it, `PDBFile.writeFile` doesn't
emit CONECT for H atoms, and downstream OpenBabel/xtb proximity bonding
may miss the bond (especially for H placed >1.9 Å from parent), causing
H atoms to drift during refinement.

### `minimize.refine_with_obminimize` / `refine_with_xtb` — universal-FF post-pass

After OpenMM AMBER minimization, optional refinement using:
- **OpenBabel obminimize** (MMFF94 → UFF auto-fallback) — fast, SMARTS-rule typing
- **xtb GFN-FF** — universal FF, GNN-derived parameters

Both use `_extract_heterogen_subsystem` to build a sub-topology of
heterogens + their full anchor residues (e.g. entire ASN, not just ND2).
Anchor atoms are frozen via:
- `OBFFConstraints.AddAtomConstraint()` for obminimize (Python API only —
  CLI obminimize has no freeze flag)
- `$fix atoms:` xcontrol block for xtb (1-based ranges)

After refinement, anchor coords are NOT spliced back (they didn't move;
keep OpenMM-minimized values). Auto-switches to heterogens-only when
topology > 5000 atoms (UFF OOMs, xtb hours-slow).

### `acpype_export.prepare_for_openmm` — glycoprotein normalization

Critical glue for the ACPYPE path. Handles:
- CYS → CYX renaming from CONECT-detected SS bonds, plus distance
  fallback (SG-SG within 2.5 Å) via `detect_ss_bonds`. Recognizes SG
  on CYS, CYX, and CYM residues.
- AMBER variant capture/restore across PDBFile normalization
  (PDBFile.writeFile converts ASH→ASP, HIE→HIS, CYX→CYS on output).
  Includes CYX and CYM (not just HID/HIE/HIP/ASH/GLH/LYN) so the
  downstream `res_templates` loop can disambiguate against CYM/CYX.
- Chain reordering: protein first, glycan after.
- Glycan TER record filtering (prevents chain breaks).
- Terminal atom removal (OXT/H2/H3 from mid-chain residues).
- HD21 addition for NLN (GLYCAM expects only HD21, not HD22).
- Peptide bond repair for stretched bonds after transplant.
- **Selective H stripping** via per-residue `_GLYCAM_KEEP_H` allowlist
  (NLN keeps `{H, HA, HB2, HB3, HD21}`, OLS keeps `{H, HA, HB2, HB3}`,
  OLT keeps `{H, HA, HB, HG21/HG22/HG23}`). Only disallowed atoms are
  stripped (HD22 on NLN, HG on OLS, HG1 on OLT, CHARMM HN/HT*). This
  preserves the canonical HD21 geometry from minimize through the
  acpype pipeline.

### `acpype_export.add_glycam_bonds(topology, ff, positions=None)` — glycan bond population

OpenMM's `PDBFile` doesn't infer bonds for non-standard residues.
This helper adds three classes of bonds for GLYCAM residues:
1. **Intra-residue bonds** — read from the FF template's bond list and
   added by atom-name match.
2. **Peptide bonds for NLN/OLS/OLT** — `prev.C → this.N` and
   `this.C → next.N`.
3. **Sugar-sugar glycosidic bonds** — distance-based detection (only
   when `positions` is provided). Anomeric C1 (or C2 for sialic) within
   2.0 Å of a linkage O2/O3/O4/O6 on another sugar. Without these,
   GLYCAM templates for linkage-position sugars (e.g. 6LB declares O6
   externally bonded) fail to match in `createSystem` ("missing 1
   externally bonded O atom"). All three of prepare/minimize/protonate
   call this with `positions=modeller.positions`.

### `transplant.py` — graft workflow

Designed for the GLYCAM-Web glycoprotein workflow: extract glycosylation
site residues from acceptor → submit to GLYCAM-Web → transplant graft
(NLN/OLS/OLT + glycan trees) back. Kabsch superposition aligns donor→
acceptor; same transform applied to graft. CONECT records remapped via
atom identity (chain, resseq, atomname) — NOT serial numbers (collide
between sources).

`--relax` runs 4-stage AMBER14+GLYCAM minimization with decreasing
restraints on protein heavy atoms (k=1000→100→10→0 kJ/mol/nm²).
`--gromacs DIR` exports GMX topology via the same `acpype_export.py`
helpers used by `dvbfixer top --acpype`.

### `glycam.convert_to_charmm(input_path, output_path)` — reverse direction

Mirror of `convert_to_glycam`. Strips GLYCAM linkage characters, looks up
the (sugar_letter, anomer) pair in `_SUGAR_LETTER_TO_CHARMM` to get the
standard PDB 3-char code (NAG/NDG/BMA/MAN/GAL/FUC/FUL/SIA/NGA/A2G/...).
Reverts NLN/OLS/OLT → ASN/SER/THR via `GLYCAM_TO_STANDARD_PROTEIN`.
Drops ROH/OME/TBT/CMET caps (CHARMM has no reducing-end caps).

Atom-name inversion uses the dict-inverted forward maps:
- `_REV_HYDROXYL_H` (universal H<n>O → HO<n>)
- `_REV_NACETYL_PDB` (NAG/NDG/NGA/A2G: C2N→C7, O2N→O7, CME→C8,
  H2N→HN2, H1M/H2M/H3M→H81/H82/H83)
- `_REV_GLYCAM_ATOM_MAP['SIA']` (sialic: C5N→C10, CME→C11, O5N→O10,
  H1M/H2M/H3M→H111/H112/H113, H5N→HN5, H3A/H3E→H31/H32, H9R/H9S→H91/H92)

Output uses 3-char PDB sugar codes that work natively with both
CHARMM-GUI (as input) and `dvbfixer top --ff charmm` (which maps them
to CHARMM RTP names via `PDB_TO_CARB`). Linkage information is
preserved via CONECT records — the residue name no longer carries it.

CLI: `dvbfixer glycam <input> --to-charmm -o <output>`.

### `prepare._preprocess_glycoprotein_input` — pre-PDBFixer input fixes

Runs before `_canonicalize_conect_records` in `run_pdbfixer`. Two
edits on the input PDB text:

1. **HETATM → ATOM rewrite** for any residue in `FORCE_ATOM_RESIDUES`
   (20 std AA + AMBER variants HID/HIE/HIP/ASH/GLH/CYX/CYM/LYN +
   GLYCAM glycoprotein residues NLN/OLS/OLT). HETATM gets treated as
   ligand by OpenMM, breaking peptide bond inference to neighbours.

2. **Spurious TER record removal** between two amino-acid residues on
   the SAME chain. A TER forces OpenMM to split the chain, breaking
   the polymer. Implementation buffers each TER and emits it only
   after seeing the next ATOM/HETATM — drops if both flanking residues
   are protein on the same chain.

Both edits are no-ops on clean inputs (returns the original path,
no temp file created).

### `prepare.find_glycosylated_atoms_with_sugar` — FF-agnostic detection

Returns `dict {(chain, resid, atom): bonded_sugar_resname}`. Two passes:

1. **CONECT-based**: scan CONECT records for bonds where one end is
   a protein anchor atom (ASN ND2, SER OG, THR OG1) and the other is
   a sugar (anomeric C of any name in `SUGAR_RESNAMES` or matching
   `is_glycam_sugar`).

2. **Distance-based fallback**: for each protein anchor atom not yet
   detected, find the nearest sugar anomeric C (C1, or C2 for sialic)
   within 2.0 Å. Catches inputs without CONECT records.

The bonded sugar name lets `rename_glycosylated_protein_residues`
branch on FF: ASN→NLN renaming fires only when `is_glycam_sugar(name)`
is True. For PDB/CHARMM sugars, ASN keeps its standard name; HD22
removal still happens via `remove_extra_glycan_hydrogens`.

### `ffutils.py` — shared GLYCAM helpers

All four GLYCAM-aware tools (`prepare`, `minimize`, `protonate`,
`top --acpype`) share these utilities:

- `is_glycam_sugar(name)` → True for GLYCAM 3-char sugar codes
  (linkage `[0-9VWUZXYTSRQPvwuzxytsr]` + sugar + anomer `[AB]`, e.g.
  UYB/4YB/VMB/0YA) or caps (ROH/OME/TBT/CMET).
- `is_glycam_residue(name)` → True for any GLYCAM sugar OR glycoprotein
  residue (NLN/OLS/OLT).
- `detect_glycam_input(topology)` → returns
  `{'glycam_proteins', 'glycam_sugars', 'pdb_sugars', 'unknown_hets'}`
  sets of (chain, res_id). Drives the FF auto-swap in `protonate` and
  the H-polish short-circuit in `prepare`.
- `fix_atom_hetatm_records(pdb_path)` → post-processes the output PDB
  to rewrite `HETATM` → `ATOM` for any residue in `FORCE_ATOM_RESIDUES`
  (= `PROTEIN_RESIDUES`, which includes the 20 std + AMBER variants
  HID/HIE/HIP/ASH/GLH/CYX/CYM/LYN + GLYCAM glycoprotein residues
  NLN/OLS/OLT). OpenMM's `PDBFile.writeFile` defaults non-standard
  residue names to HETATM; this restores them to ATOM. Idempotent.
  Called after every PDB write in the three pipeline tools.
- Constants: `GLYCAM_PROTEIN_RESIDUES = {'NLN', 'OLS', 'OLT'}`,
  `GLYCAM_CAPS = {'ROH', 'OME', 'TBT', 'CMET'}`,
  `FORCE_ATOM_RESIDUES = frozenset(PROTEIN_RESIDUES)`.

## Cross-module patterns

### CLI dispatch (`cli.py`)

Each subcommand exposes `parse_args(argv=None)` and `main(argv=None)`.
`cli.py` reads `sys.argv[1]` (command name) and dispatches `sys.argv[2:]`
to the matching module's `main()`. Entry point in `pyproject.toml`:
`dvbfixer = "dvbfixer.cli:main"`.

### CONECT record handling (`pdbutils.py`)

Three shared helpers:
- `build_serial_map(pdb_path)` → `(chain, resid, atomname) → serial`
- `remap_conect_records(input, new_serial_map)` → list of remapped CONECT lines
- `append_before_end(path, extra_lines)` → splice lines before END record

Used by `prepare`, `transplant`, `pull`, `top` for cross-stage bond
preservation. `prepare._canonicalize_conect_records` fixes malformed
CONECT spacing (5-digit serials that don't fit fixed-width columns).

### `.dat` file format (prepare/model/minimize)

JSON. Schema:
```json
{
  "description": "...",
  "total_added": 142,
  "residue_summary": {"A/GLY105": {"heavy": 4, "hydrogen": 3}},
  "added_atoms": [
    {"chain": "A", "resid": "105", "icode": "", "resname": "GLY",
     "atom": "N", "element": "N"}
  ],
  "variant_overrides": {"A:83": "HIP"}
}
```

`prepare` merges upstream `.dat` (from `model`) with its own additions
via key matching on `(chain, resid, icode, atom)`. `minimize` uses the
merged `.dat` to classify atoms into restraint tiers.

### Topology rebuild pattern (`prepare`, `minimize`)

When new atoms are added to an OpenMM topology (e.g. H atoms from RDKit),
the pattern is to build a fresh `Topology()` rather than mutating the
existing one:

1. Iterate `chain → res → atom` of the original topology
2. `new_top.addChain`, `addResidue`, `addAtom` for each
3. For each parent atom with new H to insert: `addAtom` the H, then
   `addBond(parent, h_atom)` to preserve bond connectivity
4. Re-add all `topology.bonds()` (mapped through old→new atom dict)
5. Wrap positions in `Quantity(list_of_Vec3, nanometer)`

Always include `addBond` for new H atoms — see "Heterogen H drift" gotcha
in BEST_PRACTICES.md.

## Recommended areas for future work

1. **Per-ligand antechamber for arbitrary ligands** — currently `minimize`
   falls back to strip-and-splice when SMIRNOFF can't match a residue.
   Adding `GAFFTemplateGenerator` integration (already in
   `openmmforcefields`) would unlock whole-system minimization for any
   ligand, at the cost of antechamber runtime (~30 s/ligand).

2. **xtb v6.8+ when conda-forge ships it** — fixes the `$fix` bug that
   causes ~0.1 Å drift on nominally-frozen anchors.

3. **Tests** — currently the `test/` directory contains structures used
   for manual smoke testing. A pytest suite covering the main paths
   (pure protein, GLYCAM glycoprotein, PDB-sugar glycoprotein, ligand-
   only) would catch regressions.

4. **Stand-alone modular monolith** — `top.py` (3500 lines) and
   `prepare.py` (1340 lines) carry most of the FF complexity. They're
   internally cohesive but would benefit from splitting `top.py` into
   `top/{rtp_amber.py, rtp_charmm.py, glycan.py, glycolipid.py,
   small_mol.py, output.py}` if it grows further.

5. **Docs** — README is good but long. Consider splitting per-command
   man pages or moving examples to `docs/`.
