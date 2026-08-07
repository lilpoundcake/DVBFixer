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

A working Modeller installation and license are required.
