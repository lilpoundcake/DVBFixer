# dvbfixer puppet — Backbone-Only Polyglycine Model

[← command index](index.md) · [← README](../../README.md)

Strips a PDB to a minimal backbone scaffold: removes all non-ATOM lines, removes all sidechain and hydrogen atoms (keeps only N, CA, C, O, OXT), and renames every residue to GLY. Useful for creating "puppet" models for backbone-level alignment, modeling templates, or visualization.

## Usage

```bash
# Basic usage — writes input_puppet.pdb
dvbfixer puppet input.pdb

# Keep specific residues intact (all atoms, original name)
dvbfixer puppet input.pdb --keep A:307

# Keep a range and a list
dvbfixer puppet input.pdb --keep H:286-296 --keep K:307,309

# Custom output
dvbfixer puppet input.pdb -o backbone.pdb
```

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o`, `--output` | `<input>_puppet.pdb` | Output file path |
| `--keep` | none | Keep residue(s) intact: CHAIN:NUM, CHAIN:START-END, or CHAIN:NUM1,NUM2,START-END (repeatable) |

## See also

- [Command index](index.md)
