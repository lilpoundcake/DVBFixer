# dvbfixer protonate — PROPKA3 Protonation Assignment

[← command index](index.md) · [← README](../../README.md)

Runs PROPKA3 to predict per-residue pKa values, then renames titratable residues to their correct protonation state at the target pH and adds the corresponding hydrogen atoms. Uses AMBER force field naming conventions (HID/HIE/HIP, ASH, GLH, CYM, CYX, LYN). Existing hydrogens are stripped and re-added by OpenMM based on the renamed residue templates.

**GLYCAM support**: PROPKA3 doesn't recognize GLYCAM glycoprotein residues (NLN/OLS/OLT), so they are renamed to ASN/SER/THR in a temp PDB (heterogens stripped) before PROPKA runs; pKa results are mapped back to the original input. When GLYCAM residues are detected, the FF auto-switches to AMBER14 + GLYCAM_06j-1 + tip3pfb (ff19SB has no GLYCAM templates → crash); `add_glycam_bonds(positions=...)` populates the missing bonds and `glycam-hydrogens.xml` provides the H definitions. PROPKA renames are filtered out for NLN/OLS/OLT positions (sugar-bonded sidechains, different chemistry). AMBER variant names (HID/HIE/HIP/CYX/etc.) already present in the input are preserved by scanning the raw PDB text before OpenMM normalizes them. After writing, `fix_atom_hetatm_records` rewrites HETATM lines for protein residues back to ATOM.

## Protonation Logic

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

## Usage

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

## Options

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

## See also

- [`prepare`](prepare.md) — run before `protonate` for structure repair
- [`minimize`](minimize.md) — re-runs after `protonate` to relax hydrogens
- [`rename`](rename.md) — reverse AMBER names back to canonical PDB names
