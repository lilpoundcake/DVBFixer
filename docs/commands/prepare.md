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
