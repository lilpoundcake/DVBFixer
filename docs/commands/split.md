# dvbfixer split — Empirical Chain Splitting

[← command index](index.md) · [← README](../../README.md)

Splits chains in PDB or GRO files that lack chain IDs (e.g. GROMACS MD output). Assigns unique chain IDs (A-Z, a-z, 0-9), inserts TER records, and renumbers residues per chain. GRO files are converted to PDB via MDAnalysis, preserving all residue names including protonation variants (GLUP, ASPP, etc.). Water, ions, and buffer particles (BUF/BUFF) are removed before chain detection to prevent false breaks, then optionally re-appended.

Chain breaks are detected by three criteria, applied in priority order:

1. **Residue number backward jump** — residue sequence number decreases (insertion codes like 82->82A are handled correctly and NOT treated as breaks)
2. **C->N peptide bond distance** — distance exceeds 2.5 A (any residue with backbone C/N atoms)
3. **Nearest-atom gap** — minimum distance between any atoms of consecutive residues exceeds 15 A (fallback for sugars, ligands, ions that lack peptide bonds)

**Multi-MODEL inputs** (multi-state PDBs, NMR ensembles, GROMACS trajectory exports with MODEL records) are handled as one complex sampled at multiple states: every MODEL gets the SAME chain IDs (A, B, C in every MODEL — not A B C / D E F / G H I as a naive walk would produce). The per-MODEL chain signature (atom count + residue count + first/last resname per chain) is compared across MODELs; when all match, chain IDs are reused. If MODELs differ structurally the tool falls back to independent per-MODEL chain IDs with a warning. Atom serials reset within each MODEL (standard PDB convention).

**Small-molecule threshold (`--max-chains`)** — when more than N chains are detected (default 26), only **protein** chains get chain IDs. Small-molecule chains (ions, ligands, lipids, single-residue HETATMs, glycan trees) keep a blank chain ID. A chain is classified "protein" when ≥50% of its residues are standard amino acids (incl. AMBER protonation variants HID/HIE/HIP/ASH/GLH/CYX/CYM/LYN, GLYCAM glycoprotein NLN/OLS/OLT, ACE/NME caps, MSE). Useful for structures with dozens of crystallographic ions/ligands/lipids that don't need their own letter. Raise the threshold with `--max-chains 62` to keep the original assign-all behaviour (cap is `len(CHAIN_IDS) = 62`: A-Z + a-z + 0-9).

## Usage

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

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o`, `--output` | `<input>_split.pdb` | Output file path |
| `-d`, `--distance-cutoff` | 2.5 | C->N peptide bond cutoff (angstroms) |
| `-g`, `--gap-cutoff` | 15.0 | Nearest-atom gap cutoff for non-protein residues (angstroms) |
| `--no-distance` | off | Disable all distance-based detection |
| `--no-renumber` | off | Keep original residue numbers |
| `--keep-water` | off | Keep water and ions in output (removed by default) |
| `--max-chains` | 26 | Above this many detected chains, small-molecule chains (ions, ligands, lipids, single-residue HETATMs, glycan trees) get blank chain ID; only protein chains get IDs. |
| `-v`, `--verbose` | off | Print detected chain info |

## See also

- [`renumber`](renumber.md) — SEQRES-based residue renumbering after chain splitting
- [`model`](model.md) — rebuild missing loops in the split structure
- [`zbs`](zbs.md) — full pipeline that includes splitting + renumbering + modeling
