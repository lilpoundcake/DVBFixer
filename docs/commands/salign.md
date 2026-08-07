# dvbfixer salign — Structural Alignment

[← command index](index.md) · [← README](../../README.md)

Create a structure-guided multiple alignment using the Modeller dependency.
Inputs may be complete PDB files or a single chain written as `PATH:CHAIN`.
The primary output is a Modeller PIR alignment; `--fit-dir` also retains the
superposed structures produced by SALIGN.

```bash
dvbfixer salign template1.pdb:A template2.pdb:H \
  -o templates.pir --fit-dir fitted -v
```

## Supported options

### Input and output

| Key / argument | Value | Default | Description |
|---|---|---|---|
| `template` | `PDB` or `PDB:CHAIN` | required | Two or more template structures. Add `:CHAIN` to restrict an input to one chain. |
| `-o`, `--output` | path | `structural_alignment.pir` | Output path for the Modeller PIR structural alignment. |
| `--fit-dir` | directory | none | Retain the fitted/superposed PDB structures in this directory. Without it, fitted structures are temporary. |

### Structural alignment

| Key | Value | Default | Description |
|---|---|---|---|
| `--fit-atoms` | Modeller atom selection | `CA` | Atom type or selection passed to SALIGN for structural fitting. |
| `--rms-cutoff` | ångströms | `3.5` | RMS-distance cutoff passed to SALIGN while improving the alignment. |

### Diagnostics

| Key | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` | flag | off | Enable verbose Modeller logging and SALIGN quality output. |
| `-h`, `--help` | flag | off | Print command help and exit. |

A working Modeller installation and license are required.
