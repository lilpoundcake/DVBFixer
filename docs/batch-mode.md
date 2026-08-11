# Batch mode

Batch mode runs one supported single-structure command independently for every
`.pdb` or `.ent` file in a directory (`split` also accepts `.gro`). Outputs are
written under a separate directory; `--recursive` preserves relative paths and
processing continues after individual failures unless `--fail-fast` is used.

```bash
dvbfixer zbs --input-dir structures --output-dir fixed --recursive
```

The shared keys are `--input-dir DIR`, `--output-dir DIR`, `--recursive`, and
`--fail-fast`. Supported tools are `split`, `renumber`, `model`, `pull`,
`prepare`, `minimize`, `protonate`, `rename`, `convert`, `conect`, `puppet`,
`diagnose`, and `zbs`.
