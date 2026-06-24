# `dvbfixer parametrize` — GAFF2 small molecule parametrization

[← Command index](index.md) · [← README](../../README.md)

Parametrizes an arbitrary small molecule with **GAFF2** force field (single
point of compatibility with AMBER protein force fields) and writes a
ready-to-use GROMACS topology bundle: `<name>.itp`, `<name>.gro`,
`posre_<name>.itp`. Designed for the parts of a system that don't have
RTP/template coverage in `dvbfixer top` — buffer components, drug molecules,
co-factors, organic ligands.

## Pipeline

```
input (.pdb/.mol2/.sdf)
    └─→ antechamber  (GAFF2 atom types + charges)
            └─→ parmchk2  (fill in any missing bonded params)
                    └─→ tleap  (AMBER prmtop + inpcrd)
                            └─→ ParmEd  (AMBER → GROMACS .top + .gro)
                                    └─→ dvbfixer  (split into .itp + posre.itp)
```

All AmberTools binaries ship with the bundled environment (`antechamber`,
`parmchk2`, `tleap`); no separate install needed.

## Charge methods

| Method | Flag | Speed | Quality | When to use |
|--------|------|-------|---------|-------------|
| **AM1-BCC** | `-c bcc` (default) | Seconds | Within ~5% of RESP for organic molecules. AMBER's recommended fast alternative. | Most cases. Default. |
| **RESP** | `-c resp` | Minutes (Gaussian QM, but each calculation costs as much as a small protein QM job) | Highest quality, AMBER-standard methodology | Buffer ions, conjugated systems, charged molecules where you care about the dipole. |

If unsure: **start with AM1-BCC**. Switch to RESP only when you have
empirical reason to think the BCC charges are off (compare against
published RESP charges for a similar molecule).

## Quick start (AM1-BCC)

One command, ~10 seconds:

```bash
# Neutral molecule
dvbfixer parametrize molecule.pdb -n MOL -v

# Charged: e.g. acetate (charge -1)
dvbfixer parametrize acetate.pdb -n ACET --net-charge -1 -v

# Open-shell radical: e.g. NO* (multiplicity 2)
dvbfixer parametrize radical.pdb -n RAD --net-charge 0 --multiplicity 2 -v
```

Output (in current directory):
- `ACET.itp` — `[ atomtypes ]` + `[ moleculetype ]` ready to `#include`
- `ACET.gro` — coordinates with FF atom names
- `posre_ACET.itp` — position restraints (`#ifdef POSRES`)

## RESP workflow (3 steps)

RESP requires running Gaussian, which is not bundled. You need a Gaussian
license + installation (g09 or g16).

### Step 1 — generate the Gaussian input

```bash
dvbfixer parametrize molecule.pdb -n MOL --net-charge -1 \
    -c resp --gen-gaussian
```

Writes `molecule.com`. The file is a complete Gaussian input with:

```
%mem=4GB
%nproc=4
%chk=molecule.chk
--Link1--
#HF/6-31G* SCF=tight Test Pop=MK iop(6/33=2) iop(6/42=6) opt
# iop(6/50=1)

RESP charges for molecule (q=-1, m=1) — YYYY-MM-DD

-1   1
    C   ...
    ...

molecule.gesp
```

This is the AMBER-standard RESP recipe: **HF/6-31G\*** with Merz-Kollman
ESP grid (`Pop=MK`, `iop(6/33=2)` writes ESP to log, `iop(6/42=6)`
sets grid density). Geometry is optimised first (`opt`) so charges
reflect the equilibrium structure. The trailing `molecule.gesp` writes
the ESP grid to a separate file consumed by antechamber later.

Resource controls:
- `--gaussian-mem 8GB` — for larger molecules (>50 atoms)
- `--gaussian-nproc 8` — parallel cores
- `--gaussian-method 'B3LYP/6-31G*'` — non-default method (rarely needed;
  changes the resulting charges)

### Step 2 — run Gaussian

On your QM-capable machine:

```bash
g16 < molecule.com > molecule.log
# or
g09 < molecule.com > molecule.log
```

Both `g09` and `g16` produce the same antechamber-compatible output.
Runtime depends on molecule size: a 20-atom organic typically finishes
in 10-30 minutes on 4 cores.

Verify success:
```bash
tail molecule.log    # should end with "Normal termination of Gaussian"
ls molecule.gesp     # ESP grid file should exist
```

### Step 3 — finish parametrization

Pass the `.log` back to dvbfixer:

```bash
dvbfixer parametrize molecule.pdb -n MOL --net-charge -1 \
    -c resp --gaussian-log molecule.log
```

Antechamber extracts the ESP from the log, fits RESP charges, and
emits the same `.itp` / `.gro` / `posre_*.itp` triplet as the AM1-BCC
path.

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `input` | (required) | Input structure: `.pdb`, `.mol2`, `.sdf`, or `.mol` |
| `-o`, `--output` | input stem | Output file prefix |
| `-n`, `--name` | input stem (uppercased) | Name for `[ moleculetype ]` (≤4 chars works best) |
| `-c`, `--charge-method` | `bcc` | `bcc` (AM1-BCC) or `resp` |
| `--net-charge` | `0` | Net molecular charge (integer) |
| `--multiplicity` | `1` | Spin multiplicity (1 = singlet, 2 = doublet, …) |
| `--gaussian-log` | none | Gaussian output `.log` (RESP only) |
| `--gen-gaussian` | off | Generate Gaussian input `.com` and exit (RESP only) |
| `--gaussian-mem` | `4GB` | `%mem=` directive in generated `.com` |
| `--gaussian-nproc` | `4` | `%nproc=` directive in generated `.com` |
| `--gaussian-method` | `HF/6-31G*` | QM method (don't change unless you know why) |
| `--keep-intermediate` | off | Keep `mol.mol2`, `mol.frcmod`, `mol.prmtop`, `leap.log` |
| `-v`, `--verbose` | off | Show every antechamber/tleap/parmchk2 invocation |

## Input format notes

- **PDB** — atoms inferred from records; element from column 76-78 (or atom name).
  Hydrogens must already be present (run `dvbfixer prepare` first if needed).
- **MOL2** — preferred for mixed valence and complex aromatic systems
  (bond orders are explicit).
- **SDF** — supports stereo perception and is the easiest to round-trip
  from RDKit / OpenBabel.

For PDBs without proper coordinates (e.g. straight from a SMILES + Open
Babel `--gen3d` invocation), `obminimize -ff UFF` cleanup before
parametrize improves antechamber atom-typing reliability:

```bash
obabel input.pdb -O cleaned.pdb --minimize --ff UFF
dvbfixer parametrize cleaned.pdb -n MOL
```

## Output files

For `-n MOL` and the default AM1-BCC flow, three files in the output
directory:

| File | Purpose |
|------|---------|
| `MOL.itp` | `[ defaults ]` + `[ atomtypes ]` + `[ moleculetype ]` |
| `MOL.gro` | Coordinates with FF atom names |
| `posre_MOL.itp` | Position restraints for `#ifdef POSRES` blocks |

The `.itp` includes its own `[ atomtypes ]` block, so it can be `#include`-d
**before** the protein moleculetype in a `topol.top`. Example combined
with `dvbfixer top` output:

```
#include "ffparams.itp"        ; from dvbfixer top
#include "MOL.itp"              ; from dvbfixer parametrize
#include "Protein_chain_A.itp" ; from dvbfixer top
#include "water.itp"
#include "ions.itp"

[ system ]
my_complex

[ molecules ]
Protein_chain_A   1
MOL               2     ; two copies of the ligand
SOL              12345
NA                  17
CL                  15
```

## Integration with `dvbfixer top`

`dvbfixer top` builds protein/glycan/lipid topology but can't parametrize
arbitrary organic molecules. The pattern:

```bash
# 1. Parametrize each ligand / buffer component separately
dvbfixer parametrize acetate.pdb -n ACE --net-charge -1
dvbfixer parametrize acetic_acid.pdb -n ACA
dvbfixer parametrize my_drug.pdb -n DRG --net-charge 0 -c resp \
    --gaussian-log my_drug.log

# 2. Build the protein topology
dvbfixer top complex.pdb --water opc -o gmx/

# 3. Combine: copy the small-molecule .itp files into gmx/, add the
#    #include lines and [ molecules ] entries to gmx/topol.top.
```

## Common gotchas

- **Net charge MUST match the formal charge of the molecule.** Antechamber
  will silently accept a wrong `--net-charge` and produce non-integer
  net topology charge. Always check the `Charge:` line in the summary
  output — it should be very close to your `--net-charge` value
  (within 1e-3 e).

- **Molecules with metals** (Zn, Fe, Mg coordination) aren't well handled
  by GAFF2 — use dedicated bonded models (MCPB.py) or non-bonded ion
  parameters from `dvbfixer top --ion-set`.

- **Aromatic perception** can fail on unusual heterocycles. If parmchk2
  reports missing parameters, inspect the `.frcmod` (use
  `--keep-intermediate`) and consider providing a `.mol2` input with
  explicit bond orders.

- **PDB hydrogens must be present.** Antechamber doesn't add hydrogens;
  it complains about open valences. Run `dvbfixer prepare` first if your
  input is a heavy-atom-only structure.

- **Gaussian normal termination** is REQUIRED for the RESP path.
  `tail molecule.log` should show `Normal termination of Gaussian`. If
  it failed (SCF convergence, memory, etc.), antechamber will fall back
  silently to AM1-BCC — verify with `grep "Charge Method:" mol.mol2`
  (should say `RESP` not `AM1-BCC`).

- **RESP charges are conformation-dependent.** The geometry that
  Gaussian optimises (and therefore the ESP) reflects the input
  structure. For molecules with flexible torsions, RESP averaged over
  multiple conformers (RESP-C1) is more rigorous but requires multiple
  Gaussian jobs — out of scope for this tool.

- **`-n NAME` truncation.** GROMACS allows 4-character moleculetype names;
  longer names are truncated by `gmx grompp`. Pick names ≤4 chars.

## See also

- [`top`](top.md) — full system topology; combine its output with
  parametrize's `.itp` files
- [`prepare`](prepare.md) — H repair if the small molecule input has
  no hydrogens
- [Pipelines / small molecule](../pipelines.md#small-molecule-parametrization)

## References

- Wang J, Wolf RM, Caldwell JW, Kollman PA, Case DA. *J Comput Chem*
  25, 1157 (2004) — AM1-BCC and original GAFF.
- Wang J, Wang W, Kollman PA, Case DA. *J Mol Graph Model* 25, 247
  (2006) — antechamber and GAFF2.
- Bayly CI, Cieplak P, Cornell W, Kollman PA. *J Phys Chem* 97, 10269
  (1993) — original RESP methodology.
- Cornell WD, Cieplak P, Bayly CI, Kollman PA. *J Am Chem Soc* 115,
  9620 (1993) — application of RESP to amino acid charges.
