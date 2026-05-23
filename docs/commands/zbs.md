# dvbfixer zbs — Full Pipeline

[← command index](index.md) · [← README](../../README.md)

Runs the complete preparation workflow in one command: **renumber → model → prepare → minimize → protonate → minimize**. Intermediate files are cleaned up by default — use `--keep-interim` to preserve them. Two minimize passes ensure correct protonation: the first keeps existing hydrogens (default) to get good heavy-atom positions, then protonate assigns AMBER protonation names (HIE/GLH/CYX etc.) based on PROPKA pKa predictions, then the second minimize uses `--rebuild-h` to strip and re-add hydrogens matching the correct protonation state (e.g. HE2 for GLH). The final output has AMBER protonation names — use `dvbfixer rename` if you need canonical PDB names. The `.dat` file flows from model (gap atoms) through prepare (merged with PDBFixer additions) to minimize (selective restraints). Water is removed by default. Each step can be skipped individually.

## Usage

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

## Options

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
| `--strip-heterogens` | off (default: keep) | Strip heterogens during prepare/minimize — protein-only pipeline |
| `--mutate` | none | Mutate a residue during prepare: CHAIN:RESNUM:NEW_AA (repeatable) |
| `--skip-minimize` | off | Skip minimize step |
| `--no-solvent` | off | Minimize in vacuum |
| `--rebuild-h` | off | Strip and re-add hydrogens via OpenMM during minimization (default: keep existing) |
| `--restraint-k` | 100.0 | Restraint force constant |
| `--max-iter` | 1000 | Max minimization iterations per phase |
| `--platform` | auto | OpenMM platform |
| `--skip-protonate` | off | Skip protonate step |
| `--no-hydrogens` | off | Only rename residues, skip hydrogen addition |
| `--keep-water` | off | Keep water molecules (removed by default) |
| `--keep-interim` | off | Keep all intermediate files (default: only final output) |
| `-v`, `--verbose` | off | Print detailed progress for all steps |

## See also

- [`model`](model.md), [`prepare`](prepare.md), [`minimize`](minimize.md), [`protonate`](protonate.md) — the individual steps
- [Pipelines](../pipelines.md) — manual step-by-step recipes and alternatives
