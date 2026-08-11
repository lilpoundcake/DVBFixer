# dvbfixer rename — Canonicalize Residue Names

[← command index](index.md) · [← README](../../README.md)

Renames non-canonical residues to standard PDB three-letter codes. Handles AMBER protonation names (HIE/HID/HIP→HIS, ASH→ASP, GLH→GLU, CYX/CYM→CYS, LYN→LYS), CHARMM names (HSD/HSE/HSP→HIS), and selenomethionine (MSE→MET). Text-based — does not modify coordinates or atoms.

## Usage

```bash
# Basic usage — writes input_renamed.pdb
dvbfixer rename input.pdb

# Verbose output
dvbfixer rename input.pdb -v
```

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o`, `--output` | `<input>_renamed.pdb` | Output file path |
| `-v`, `--verbose` | off | Print each rename |

Also available as `--rename` flag on `renumber`, `prepare`, `minimize`, and `pull` tools to canonicalize input before processing.

## See also

- [`renumber`](renumber.md) — `--rename` integrates rename into the renumber step
- [`prepare`](prepare.md) — `--rename` integrates rename into the prepare step
- [`convert`](convert.md) — bidirectional PDB↔GLYCAM sugar nomenclature (formerly `glycam`)

## How it works
Text-based rename of non-canonical residue names to standard PDB names. Converts AMBER protonation (HIE/HID/HIP→HIS, ASH→ASP, GLH→GLU, CYX/CYM→CYS, LYN→LYS), CHARMM (HSD/HSE/HSP→HIS), and MSE→MET. Also available as `--rename` flag on `renumber`, `prepare`, `minimize`, and `pull` tools. Uses `CANONICAL_MAP` dict in `rename.py`.

## Batch mode

`rename` can canonicalize a directory without combining structures:
`dvbfixer rename --input-dir structures --output-dir canonical`.
See [Batch mode](../batch-mode.md) for shared keys.
