# dvbfixer prepare — Structure Fixing with PDBFixer

[← command index](index.md) · [← README](../../README.md)

Adds missing residues, missing heavy atoms, and hydrogens using PDBFixer. **Heterogens (sugars, ligands) are kept and protonated by default** — H is added to them via AMBER14 + GLYCAM_06j-1 (BioLuminate-style), so a crystal structure with bare glycans/ligands becomes fully protonated in one call. Pass `--strip-heterogens` to remove them (protein-only mode), or `--no-heterogen-h` to keep heterogens but skip the H addition. For arbitrary unknown ligands that need real force-field parameters (not just H placement), run `dvbfixer minimize --parametrize-ligands` after prepare. User-specified protonation variants (HIE/HID/HIP/ASH/GLH/CYX from input PDB or `--mutate`) are passed as explicit OpenMM variants for correct hydrogen placement. Writes a `.dat` file recording which atoms were added (including variant overrides for downstream minimize).

**GLYCAM-named input** (NLN/OLS/OLT + GLYCAM sugars like UYB/4YB/VMB/0YA) is auto-detected and handled natively — `add_glycam_bonds(positions=...)` is called before addHydrogens to populate intra-residue, peptide, and sugar-sugar glycosidic bonds (required for template matching on NLN's protein neighbour). The RDKit/OpenBabel H-polish passes are skipped when all heterogens are GLYCAM-named (already protonated by GLYCAM templates). After writing the output, `fix_atom_hetatm_records` rewrites any HETATM lines for protein residues (HID/HIE/HIP/ASH/GLH/CYX/CYM/LYN/NLN/OLS/OLT) back to ATOM records.

**Input preprocessing** — `_preprocess_glycoprotein_input` runs first to fix two common upstream-tool issues that break OpenMM topology parsing: (1) HETATM lines for protein/GLYCAM glycoprotein residues are rewritten to ATOM (HETATM gets treated as ligand → no peptide bond inferred to neighbours → "TYR missing externally bonded C atom" template errors), and (2) spurious TER records between two amino-acid residues on the same chain are dropped (a TER forces a new chain in OpenMM, breaking the polymer). Both edits are no-ops on clean inputs.

**Glycosylation detection is FF-agnostic** — `find_glycosylated_atoms_with_sugar` uses CONECT records AND a distance-based fallback (ASN ND2 / SER OG / THR OG1 within 2.0 Å of a sugar anomeric C). Inputs with PDB sugars (NAG/NDG/BMA/...) and inputs with CHARMM-GUI 4-char sugars (BGLC/BMAN/AMAN/BGLCNA/...) are both recognized. The ASN→NLN rename fires ONLY when the bonded sugar is GLYCAM-named — for PDB/CHARMM sugars, ASN stays as ASN. The extra HD22 on glycosylated ND2 is removed in all three FF conventions (consistent across CHARMM, AMBER, GLYCAM).

**Mutations — substitutions and deletions** via `--mutate`. Two forms are supported:

- `CHAIN:RESNUM:NEW_AA` — substitute (e.g. `A:39:ALA`, `A:83:HIP` for a protonation variant)
- `CHAIN:RESNUM:del` — DELETE the residue entirely (e.g. `H:446:del`, `H:100A:del`)

Both forms run a **dependency cleanup** at the raw-PDB-text level BEFORE PDBFixer sees the structure, so the new residue (or empty slot) is internally consistent. Edge cases handled automatically:

- **Glycan removal on glycosylation-site change** — if the substituted/deleted residue's sidechain anchor (ASN ND2, SER OG, THR OG1, TYR OH, CYS SG, LYS NZ, ARG NH2) carries a HETATM tree via CONECT, the entire glycan tree is walked and removed. Applies to both deletions (`--mutate H:297:del`) and substitutions to a residue that can't carry the glycan (`--mutate H:297:ALA`). Protonation-variant renames where the parent residue is unchanged (`--mutate A:39:CYX`, `--mutate A:83:HIP`) do NOT trigger cleanup.
- **Disulfide partner repair** — substituting/deleting a CYS/CYX automatically reduces its SS partner: CYX→CYS rename + drop existing HG (so `addHydrogens` regenerates it). SSBOND records referencing the affected residue are dropped. `--mutate A:39:CYX` keeps the bridge (no parent change).
- **LINK partners** — LINK records mentioning an affected residue are dropped and a warning is printed. The partner is left as-is — inspect manually if it needs repair.
- **Terminal vs internal (deletion only)** — terminal deletions need no peptide-bond reconnection; internal deletions are bridged by the downstream `minimize` step. If the post-deletion C(i-1)→N(i+1) distance exceeds 5 Å a warning is printed (consider running `dvbfixer pull` or `dvbfixer model`).
- **Insertion codes** — `H:100A:del` selects residue 100A specifically.
- **Multiple deletions** — consecutive deletions are treated as one contiguous gap.
- **Mixed with substitutions** — deletion and substitution of the same residue is rejected at parse time.

The `.dat` file gains a `removed_residues` field with per-mutation metadata (resname, gap type which is `substitution` for substitution cleanups, gap distance, linked glycan residues, disulfide partner repaired, `substituted_to` field for substitutions).

## Usage

```bash
# Basic usage — writes input_prepared.pdb and input_prepared.dat
dvbfixer prepare input.pdb -v

# Custom output
dvbfixer prepare input.pdb -o fixed.pdb --dat fixed.dat

# Keep crystallographic waters (heterogens already kept by default)
dvbfixer prepare input.pdb --keep-water

# Strip heterogens (protein-only mode)
dvbfixer prepare input.pdb --strip-heterogens

# Apply point mutations
dvbfixer prepare input.pdb --mutate A:39:ALA --mutate B:100:GLY -v

# Delete residues (handles attached glycans + SS partners automatically)
dvbfixer prepare input.pdb --mutate H:446:del --mutate H:447:del -v

# Knock out a glycosylation site — attached glycan tree is removed too
dvbfixer prepare glycoprotein.pdb --mutate H:297:ALA -v

# Break a disulfide — partner is auto-reduced (CYX→CYS + HG regenerated)
dvbfixer prepare antibody.pdb --mutate H:22:ALA -v

# Mix substitution and deletion
dvbfixer prepare input.pdb --mutate A:39:ALA --mutate H:446:del -v
```

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o`, `--output` | `<input>_prepared.pdb` | Output PDB file |
| `--dat` | `<output>.dat` | Restraint data file path |
| `--ph` | 7.0 | pH for hydrogen addition |
| `--keep-water` | off | Keep crystallographic waters |
| `--strip-heterogens` | off (default: keep) | Remove heterogens (sugars, ligands, ions) before processing — protein-only mode |
| `--no-heterogen-h` | off | Keep heterogens but skip H addition |
| `--ff` | `auto` | Force field for the heterogen-H addition step. Accepts a short name (`auto`, `amber`, `amber+glycam`, `charmm`, …) or explicit OpenMM XML paths. Only consulted when heterogen-H addition actually runs. See [force-fields.md](../force-fields.md). |
| `--mutate` | none | Mutate a residue: `CHAIN:RESNUM:NEW_AA` (substitution) or `CHAIN:RESNUM:del` (deletion). Insertion codes supported (`H:100A:del`). Repeatable. |
| `--rename` | off | Rename non-canonical residues (AMBER/CHARMM) to standard names before processing |
| `-v`, `--verbose` | off | Print detailed progress |

## The .dat File

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

When `--mutate ...:del` is used, the `.dat` also contains a `removed_residues` list:

```json
"removed_residues": [
  {
    "chain": "C", "resid": "64", "icode": "", "resname": "ASN",
    "removed_atoms": 8, "gap_type": "internal",
    "prev_residue": "C/ASN/63", "next_residue": "C/ASP/65",
    "gap_distance_A": 3.5,
    "linked_glycan_residues": [{"chain":"C","resid":"206","resname":"NAG"}],
    "disulfide_partner_repaired": null
  }
]
```

Downstream `minimize` reads this so it knows which residues are intentionally absent (and doesn't try to match against an upstream `.dat` that still listed them).

## See also

- [`model`](model.md) — feeds its `.dat` into `prepare`
- [`minimize`](minimize.md) — consumes the merged `.dat`
- [`protonate`](protonate.md) — apply pH-based protonation after prepare
- [BEST_PRACTICES.md](../../BEST_PRACTICES.md) — glycoprotein preparation workflow

## How it works

### Who places hydrogens?

The pipeline splits the work between PDBFixer and OpenMM:

- **PDBFixer** does all *heavy*-atom repair — missing residues (via
  SEQRES gap filling), missing atoms in existing residues, terminal
  atoms (OXT), heterogen removal, mapping non-standard residues to
  their standard replacements.
- **OpenMM's `Modeller.addHydrogens(forcefield, pH=..., variants=[...])`**
  places every hydrogen.

`prepare` deliberately **does not** call PDBFixer's
`addMissingHydrogens(pH)`. That method uses PDBFixer's own
`_describeVariant` internal, which only recognises standard PDB names
(HIS / ASP / GLU / CYS / LYS). It silently ignores AMBER variant
labels (HIE / HID / HIP / ASH / GLH / CYX / CYM / LYN) coming from
the input PDB, from `--mutate`, or from downstream PROPKA / Reduce
picks. Only `Modeller.addHydrogens` takes an explicit `variants=`
list, so that's the only entry point that respects user intent about
protonation state.

The trade-off: it's OpenMM's `addHydrogens` — not PDBFixer's — that
has the CSER-template misplacement bug (0.4.1 / 0.4.2 fixed by
`dvbfixer.ffutils.geometry.repair_misplaced_hydrogens`, run
immediately after every `addHydrogens` call). Since `prepare`,
`protonate`, and `minimize` all use OpenMM's `addHydrogens` for the
same reason, the post-check runs after each one.

### Original overview
Runs PDBFixer to add missing residues, heavy atoms, and hydrogens. **Default: keep heterogens AND add H to them** (BioLuminate-style) via `create_forcefield_with_openff(['amber14-all.xml', 'amber14/GLYCAM_06j-1.xml'], topology)` + `Modeller.loadHydrogenDefinitions('glycam-hydrogens.xml')` + `add_glycam_bonds(positions=...)` + `addHydrogens(ff, ...)`. The `add_glycam_bonds` call BEFORE addHydrogens is critical — it populates intra-residue bonds for NLN/OLS/OLT plus peptide bonds (NLN.C → next.N) plus distance-based sugar-sugar glycosidic bonds (anomeric C1/C2 within 2.0 Å of linkage O2/O3/O4/O6). Without these, template matching on the residue ADJACENT to NLN fails ("TYR missing 1 externally bonded C atom"). Wrapped in try/except → on failure (e.g. unknown ligand without SMILES), falls back to protein-only `addHydrogens(pH=ph)` call. Strips H from all non-solvent/non-ion residues before the call so they regenerate with correct GLYCAM atom names. **GLYCAM short-circuit**: when every heterogen is GLYCAM-named (`is_glycam_residue` returns True for all), the post-addHydrogens RDKit/OpenBabel H-polish passes are SKIPPED — they would strip atoms with GLYCAM-specific names (C2N/O2N/CME) that don't match OpenBabel's valence rules. CLI: `--strip-heterogens` (opt-in protein-only mode), `--no-heterogen-h` (keep heterogens but skip H addition). Bypasses PDBFixer's `addMissingHydrogens` (which ignores topology renames) — calls `Modeller.addHydrogens(variants=...)` directly with explicit variants list. User variant overrides (from `--mutate` or input PDB HIE/ASH/GLH) passed as OpenMM variants for correct H placement; standard HIS uses `None` (pH auto-detect). Output PDB written with standard names via `PDBFile.writeFile`, then text-based post-processing restores user variant names, then `fix_atom_hetatm_records` rewrites any remaining HETATM lines for protein/GLYCAM residues to ATOM. NLN/OLS/OLT are filtered OUT of PDBFixer's nonstandardResidues print BEFORE the message emits (avoids misleading "NLN -> LEU" warning that never actually happens because they're also filtered before `replaceNonstandardResidues`). Writes a `.dat` file (JSON) recording added atoms + `variant_overrides` for downstream minimize. Merges upstream `.dat` (from model step) if present. Detects glycosylated residues from CONECT records AND distance fallback and removes extra hydrogens (e.g. ASN HD22 on glycosylated ND2). Supports `--mutate CHAIN:RESNUM:NEW_AA` for point mutations including protonation variants (e.g. `--mutate A:83:HIP`). `--rename` canonicalizes non-standard residue names before processing.

**Mutation cleanup** (`--mutate CHAIN:RESNUM:del` AND `--mutate CHAIN:RESNUM:NEW_AA` when the new AA changes the parent residue): `parse_mutations` returns a third value `deletions`; `apply_deletions_to_pdb_text` runs as a new step BEFORE `_preprocess_glycoprotein_input` and operates on raw PDB text. It accepts both deletion and substitution-cleanup targets. For each target it: (1) for deletions only, marks the residue's atoms for removal; for substitution-cleanup, the residue itself stays — only its dependent atoms/records are cleaned up; (2) **glycan walk**: BFS through `CONECT` from the residue's sidechain anchor (`SIDECHAIN_ANCHORS`: ASN ND2, GLN NE2, SER OG, THR OG1, TYR OH, HYP OD1, CYS/CYX/CYM SG, LYS NZ, ARG NH2), crossing only into HETATM atoms of other residues — every reachable HETATM residue is added to the removal set, so an N-linked NAG-NAG-BMA-... tree disappears with the ASN. Applies equally to substitution (e.g. `--mutate H:297:ALA` removes the glycan because ALA can't carry it) and deletion; (3) **disulfide partner**: if the affected residue is CYS/CYX, find the partner via SSBOND records (or fallback CONECT SG-SG); the partner is renamed `CYX→CYS` and its `HG` is dropped so `addHydrogens` regenerates it; the SSBOND record is dropped. Applies equally to substitution (e.g. `--mutate H:22:ALA`) and deletion; (4) **LINK partners**: any LINK record naming the affected residue is dropped with a warning, the partner is left as-is; (5) **terminal vs internal classification (deletions only)**: by file-order position in the chain (skipping other deletions); for internal gaps, the `prev.C → next.N` distance is computed from the input coordinates and a warning emitted if > 5 Å.

Substitution-cleanup is filtered: when the user-supplied NEW_AA's standard parent name matches the old residue name (e.g. `--mutate A:39:CYX` where the input has CYS — CYX is just a protonation-variant rename, parent is still CYS), no cleanup is run. This preserves SS bonds when the user is just annotating an existing disulfide-bonded cysteine with the CYX name. Same logic applies to ASN→variant renames (none defined in AMBER) etc.

The cleaned PDB is written with: atom lines filtered (deletions only), CONECT lines rewritten without removed serials, dropped SSBOND/LINK lines elided, CYX→CYS partners renamed in-place with HG dropped. After PDBFixer loads the cleaned file, `fixer.missingResidues` entries that correspond to **deleted** positions are scrubbed between `findMissingResidues()` and `findMissingAtoms()` — otherwise PDBFixer's SEQRES-driven gap filler would re-add the residue we just removed. The scrub handles chains where multiple OpenMM `Chain` objects share a chain ID (e.g. a HETATM block at the end of the file gets its own `Chain` object even though it has the same letter as a protein chain). The `.dat` file gains a `removed_residues` field with per-mutation metadata: chain, resid, icode, resname, removed_atoms count (0 for substitutions), gap_type (internal/terminal_N/terminal_C/whole_chain/`substitution`), gap_distance_A, prev_residue/next_residue strings (deletion-only), linked_glycan_residues list, disulfide_partner_repaired, and a `substituted_to` field for substitutions. Substitution and deletion of the same residue is rejected at parse time. Substitution mutations cannot use insertion codes in RESNUM (PDBFixer's `applyMutations` doesn't address iCode residues); deletions can.

**Input preprocessing** (`_preprocess_glycoprotein_input`, runs before `_canonicalize_conect_records`): fixes two common upstream-tool issues. (1) `HETATM` lines for residues in `FORCE_ATOM_RESIDUES` (20 std AA + AMBER variants + NLN/OLS/OLT) are rewritten to `ATOM  ` — HETATM gets treated as ligand by OpenMM, breaking peptide bond inference to neighbours and producing "TYR missing externally bonded C atom" template errors. (2) Spurious TER records between two amino-acid residues on the SAME chain are dropped — a TER forces OpenMM to split the chain, breaking the polymer. Both edits are no-ops on clean inputs (returns the original path).

**FF-agnostic glycosylation detection** (`find_glycosylated_atoms_with_sugar`): returns `{(chain, resid, atom): bonded_sugar_resname}`. Uses CONECT records for protein-sugar bonds AND distance fallback (ASN ND2 / SER OG / THR OG1 within 2.0 Å of a sugar anomeric C; C2 for sialic) — catches glycosylation sites that have no CONECT record (common in CHARMM-GUI output and crystal PDBs). Sugar set includes PDB 3-char names (NAG/NDG/BMA/MAN/GAL/FUC/FUL/SIA/NGA/A2G/...), CHARMM-GUI 4-char names (BGLC/BMAN/AMAN/BGAL/BGLCNA/...), and GLYCAM 3-char codes (via `is_glycam_sugar`). The ASN→NLN rename in `rename_glycosylated_protein_residues` fires ONLY when the bonded sugar is GLYCAM-named (NLN is a GLYCAM-specific name). For PDB/CHARMM sugars, ASN/SER/THR stay with standard names; HD22/HG/HG1 removal still happens via `remove_extra_glycan_hydrogens` (consistent behavior across all three FFs).
