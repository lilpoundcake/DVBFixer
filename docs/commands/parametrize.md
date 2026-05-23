# dvbfixer parametrize — GAFF2 Small Molecule Parametrization

[← command index](index.md) · [← README](../../README.md)

Parametrizes small molecules with GAFF2 force field and AM1-BCC or RESP charges for GROMACS MD. Wraps the AmberTools pipeline: antechamber → parmchk2 → tleap → ParmEd. Output: standalone `.itp` + `.gro` + `posre.itp`.

## Usage

```bash
# AM1-BCC (default, fast)
dvbfixer parametrize molecule.pdb -n MOL -v

# Acetate with charge -1
dvbfixer parametrize acetate.pdb -n ACET --net-charge -1

# RESP (requires Gaussian log)
dvbfixer parametrize molecule.pdb -n MOL -c resp --gaussian-log molecule.log

# Generate Gaussian input for RESP
dvbfixer parametrize molecule.pdb -n MOL -c resp --gen-gaussian
```

## See also

- [`top`](top.md) — full system topology (proteins, glycans, ions); use `parametrize` for the small-molecule pieces, then `#include` the `.itp` from `topol.top`
- [`prepare`](prepare.md) — apply before `parametrize` if the small molecule needs H repair
