# dvbfixer pull — Bond Pulling

[← command index](index.md) · [← README](../../README.md)

Pulls atoms together to form bonds (disulfide bridges, glycosidic bonds) using OpenMM partial minimization. Atoms within a configurable radius of bond endpoints move freely; the rest are frozen via mass=0. Auto-removes conflicting hydrogens (CYS HG for disulfides, ASN HD22 for glycosidic bonds). Validates bond geometry before and after pulling.

## Usage

```bash
# Form a disulfide bond between CYS residues
dvbfixer pull input.pdb --bond A:22:SG:A:96:SG -v

# Multiple bonds
dvbfixer pull input.pdb --bond A:22:SG:A:96:SG --bond A:300:ND2:A:1301:C1 -v

# Custom radius and output
dvbfixer pull input.pdb --bond A:22:SG:A:96:SG --radius 8.0 -o output.pdb
```

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `--bond` | (required) | Bond specification: CHAIN1:RES1:ATOM1:CHAIN2:RES2:ATOM2 (repeatable) |
| `-o`, `--output` | `<input>_pull.pdb` | Output file path |
| `--radius` | 6.0 | Radius around bond endpoints for free atoms (angstroms) |
| `--target` | auto | Target bond distance (angstroms, auto-detected from element pair) |
| `--rename` | off | Rename non-canonical residues before processing |
| `--ff` | `auto` | Force field for the OpenMM partial minimization. Accepts a short name (`auto`, `amber`, `amber+glycam`, `charmm`, …) or explicit OpenMM XML paths. See [force-fields.md](../force-fields.md). |
| `-v`, `--verbose` | off | Print detailed progress |

## See also

- [`minimize`](minimize.md) — full-system minimization after pull
- [`transplant`](transplant.md) — alternative for assembling glycoproteins from GLYCAM-Web output

## How it works
Pulls atoms together to form bonds (disulfide bridges, glycosidic bonds) using OpenMM partial minimization. Supports multiple `--bond` specifications. Protein-protein bonds use `CustomBondForce`; protein-HETATM bonds use `CustomExternalForce` toward the fixed HETATM position. Atoms within `--radius` of bond endpoints are free to move (mass=0 freezing for the rest). Auto-removes conflicting hydrogens (CYS HG for disulfides, ASN HD22 for glycosidic bonds). Pre-pull validation checks valence (bond count vs `MAX_BONDS`) and bond type reasonableness for the pulling residues. Post-pull validation checks convergence (distance vs target), bond length range, and steric clashes within pulling residues.
