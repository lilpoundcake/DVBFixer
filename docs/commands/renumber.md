# dvbfixer renumber — SEQRES-Based Residue Renumbering

[← command index](index.md) · [← README](../../README.md)

Renumbers residues by aligning ATOM records to the SEQRES section via subsequence matching. Removes insertion codes (e.g. Kabat/Chothia antibody numbering 100A-J -> sequential 105-114) while preserving correct gap positions for missing residues.

Updates **all** PDB sections that reference residue numbers:
- ATOM, HETATM, TER
- HELIX, SHEET, SSBOND, LINK, CISPEP
- HET, DBREF, SEQADV
- CONECT (atom serial remapping)
- REMARK 465, 500, 610

## Usage

```bash
# Basic usage — writes input_renum.pdb
dvbfixer renumber input.pdb

# Verbose output showing alignment details and gaps
dvbfixer renumber input.pdb -v

# Custom output
dvbfixer renumber input.pdb -o renumbered.pdb
```

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o`, `--output` | `<input>_renum.pdb` | Output file path |
| `--keep-water` | off | Keep water molecules (HOH, WAT, TIP3, SOL) — removed by default |
| `--rename` | off | Rename non-canonical residues (AMBER/CHARMM) to standard names before processing |
| `-v`, `--verbose` | off | Print alignment details and gap positions |

## How It Works

1. Parses SEQRES records to get the full sequence per chain
2. Extracts unique (resSeq, iCode, resname) tuples from ATOM records
3. Aligns ATOM residues to SEQRES as a subsequence — each ATOM residue name is matched to the next occurrence in SEQRES
4. Assigns new sequential numbering based on SEQRES position (position 1 = first SEQRES residue)
5. Non-SEQRES residues (waters, ligands) are numbered sequentially after the last SEQRES position
6. Chains without SEQRES entries are renumbered sequentially from 1
7. All PDB sections are updated with the new numbering

## See also

- [`split`](split.md) — apply before `renumber` when the input has no chain IDs
- [`model`](model.md) — rebuild missing loops (consumes renumbered output)
- [`rename`](rename.md) — text-only rename of non-canonical residues
