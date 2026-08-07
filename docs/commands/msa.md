# dvbfixer msa — Multiple Sequence Alignment

[← command index](index.md) · [← README](../../README.md)

Align two or more protein sequences through a consistent interface. MAFFT is
the default when several engines are installed; MUSCLE 5 and Clustal Omega can
be selected explicitly. Output is validated to ensure that row identifiers and
ungapped sequences were not changed by the external program.

```bash
dvbfixer msa sequences.fasta -o alignment.fasta
dvbfixer msa sequences.fasta --engine muscle --format pir -o alignment.pir
dvbfixer msa target.fasta --template template1.pdb:A --template template2.pdb:H
dvbfixer msa ignored.fasta --list-engines
```

The executable is resolved from `PATH` as `mafft`, `muscle`, or `clustalo`.
Install the full environment to obtain all three engines, or install only the
engine required by a headless workflow. MUSCLE must be version 5; MUSCLE 3 has
a different command-line interface. Run `dvbfixer msa --list-engines` to see
the executable paths detected by the CLI. Full conda, Homebrew, Linux, Windows,
and standalone-binary instructions are in the
[installation guide](../installation.md#multiple-sequence-alignment-executables).
