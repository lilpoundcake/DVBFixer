# dvbfixer prepare — Structure Fixing with PDBFixer

[← command index](index.md) · [← README](../../README.md)

Adds missing residues, missing heavy atoms, and hydrogens using PDBFixer. **Heterogens (sugars, ligands) are kept and protonated by default** — H is added to them via AMBER14 + GLYCAM_06j-1 + SMIRNOFF (BioLuminate-style), so a crystal structure with bare glycans/ligands becomes fully protonated in one call. Pass `--strip-heterogens` to remove them (protein-only mode), or `--no-heterogen-h` to keep heterogens but skip the H addition. User-specified protonation variants (HIE/HID/HIP/ASH/GLH/CYX from input PDB or `--mutate`) are passed as explicit OpenMM variants for correct hydrogen placement. Writes a `.dat` file recording which atoms were added (including variant overrides for downstream minimize).

**GLYCAM-named input** (NLN/OLS/OLT + GLYCAM sugars like UYB/4YB/VMB/0YA) is auto-detected and handled natively — `add_glycam_bonds(positions=...)` is called before addHydrogens to populate intra-residue, peptide, and sugar-sugar glycosidic bonds (required for template matching on NLN's protein neighbour). The RDKit/OpenBabel H-polish passes are skipped when all heterogens are GLYCAM-named (already protonated by GLYCAM templates). After writing the output, `fix_atom_hetatm_records` rewrites any HETATM lines for protein residues (HID/HIE/HIP/ASH/GLH/CYX/CYM/LYN/NLN/OLS/OLT) back to ATOM records.

**Input preprocessing** — `_preprocess_glycoprotein_input` runs first to fix two common upstream-tool issues that break OpenMM topology parsing: (1) HETATM lines for protein/GLYCAM glycoprotein residues are rewritten to ATOM (HETATM gets treated as ligand → no peptide bond inferred to neighbours → "TYR missing externally bonded C atom" template errors), and (2) spurious TER records between two amino-acid residues on the same chain are dropped (a TER forces a new chain in OpenMM, breaking the polymer). Both edits are no-ops on clean inputs.

**Glycosylation detection is FF-agnostic** — `find_glycosylated_atoms_with_sugar` uses CONECT records AND a distance-based fallback (ASN ND2 / SER OG / THR OG1 within 2.0 Å of a sugar anomeric C). Inputs with PDB sugars (NAG/NDG/BMA/...) and inputs with CHARMM-GUI 4-char sugars (BGLC/BMAN/AMAN/BGLCNA/...) are both recognized. The ASN→NLN rename fires ONLY when the bonded sugar is GLYCAM-named — for PDB/CHARMM sugars, ASN stays as ASN. The extra HD22 on glycosylated ND2 is removed in all three FF conventions (consistent across CHARMM, AMBER, GLYCAM).

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
| `--mutate` | none | Mutate a residue: CHAIN:RESNUM:NEW_AA (can be used multiple times) |
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

## See also

- [`model`](model.md) — feeds its `.dat` into `prepare`
- [`minimize`](minimize.md) — consumes the merged `.dat`
- [`protonate`](protonate.md) — apply pH-based protonation after prepare
- [BEST_PRACTICES.md](../../BEST_PRACTICES.md) — glycoprotein preparation workflow
