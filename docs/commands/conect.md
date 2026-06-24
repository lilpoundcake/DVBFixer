# `dvbfixer conect`

[← Command index](index.md)

Infer missing `CONECT` records for a PDB file and write a copy with the
merged bond block. The same inference engine runs **automatically** inside
`prepare`, `top`, `minimize`, `transplant`, and `convert`, so users on a
clean pipeline rarely need to invoke this subcommand directly. Use it
when you want CONECT as a discrete preprocessing step, or to inspect
what bonds dvbfixer perceives in a given input.

## Why

Many dvbfixer tools (`prepare`'s glycosylation detection, `top`'s sugar
tree builder, `transplant`'s CONECT remapping, `minimize`'s disulfide
detection) need bond connectivity that lives in `CONECT` records. Inputs
from RCSB, GROMACS, EM depositions, and GLYCAM-Web often ship without
CONECT — and `OpenMM PDBFile.writeFile` strips CONECT silently between
pipeline stages. Auto-inference makes those flows work end-to-end.

## Algorithm

Hybrid — OpenBabel `ConnectTheDots` (Blue Obelisk standard cutoff
`r < r_cov1 + r_cov2 + 0.45 Å`) for general bond perception, then three
domain-specific overrides that match the patterns dvbfixer downstream
tools consume:

| Override | Detects | Cutoff |
|---|---|---|
| **SS** | `CYS`/`CYX`/`CYM` SG-SG | 2.5 Å |
| **Glycosidic** | Sugar anomeric C (C1, or C2 for sialic) → neighbouring sugar O2/O3/O4/O6 | 2.0 Å |
| **Glycosylation** | Protein ND2/OG/OG1 (ASN/SER/THR/NLN/OLS/OLT) → sugar anomeric C | 2.5 Å |

The output is the **union** of existing CONECT (if any) and inferred
bonds, deduplicated. Bonds where both atoms are in the same standard
amino-acid residue are dropped from the inferred set — FF templates own
that chemistry and emitting them just bloats the file.

If OpenBabel is unavailable, falls back to element-aware distance
cutoffs (KDTree-accelerated for systems ≥2000 atoms).

## Usage

```bash
# Default: write <input>_conect.pdb alongside the input
dvbfixer conect input.pdb

# Explicit output path
dvbfixer conect input.pdb -o input_with_bonds.pdb

# Verbose: print counts of (existing, inferred, total)
dvbfixer conect input.pdb -v

# In-place overwrite (rare — guarded by --force)
dvbfixer conect input.pdb -o input.pdb --force
```

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o`, `--output` | `<input>_conect.pdb` | Output PDB path |
| `--force` | off | Allow `--output` to equal input (in-place overwrite) |
| `--include-protein-backbone` | off | Also emit CONECT for standard amino-acid bonds. FF templates own those chemistry, but cyclic peptides built from standard AAs need the N-to-C closure as CONECT to be detected downstream. |
| `-v`, `--verbose` | off | Print bond counts |

## Auto-inference in other tools

The `prepare`, `top`, `minimize`, `transplant`, and `convert` tools
materialise an inferred-CONECT temp PDB at the start of `main()` and
read from that copy. **The user's input file is never modified.**

To opt out (e.g. for debugging what the input actually declared), pass
`--no-infer-conect` to any of those tools:

```bash
# Default: auto-infer runs
dvbfixer top input.pdb --acpype -o gmx/

# Disable: strict CONECT-only behaviour
dvbfixer top input.pdb --acpype -o gmx/ --no-infer-conect
```

## What's NOT inferred

- Standard amino-acid backbone/sidechain bonds — owned by FF templates.
- Bond orders, aromaticity, stereo — `CONECT` is connectivity only.
- LINK records — dvbfixer downstream uses CONECT exclusively.

## See also

- [`prepare`](prepare.md) — uses inferred glycosylation sites
- [`top`](top.md) — uses inferred sugar tree + SS bonds
- [`transplant`](transplant.md) — uses inferred bonds for CONECT remapping
- [Known issues](../known-issues.md)

## References

- Open Babel `ConnectTheDots`: O'Boyle et al., *J. Cheminformatics* 3, 33 (2011).
- Cordero covalent radii: *Dalton Trans.* DOI 10.1039/B801115J (2008).
- wwPDB CONECT format spec: https://www.wwpdb.org/documentation/file-format-content/format33/sect10.html
