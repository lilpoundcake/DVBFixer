# `dvbfixer minimize`

[User guide](../README.md) · [Command index](index.md) · [CLI reference](../reference/minimize.md)

Relaxes coordinates with OpenMM and selective restraints. Use after
[prepare](prepare.md); minimization does not replace missing-sequence modeling,
solvent equilibration, or assessment of experimental correctness.

## Quick start

```bash
dvbfixer minimize prepared.pdb --dat prepared.dat --no-solvent -o minimized.pdb
dvbfixer diagnose minimized.pdb
```

**Inputs:** a prepared structure and its matching repair sidecar.
**Result:** `minimized.pdb`. The example uses vacuum minimization for a shorter
run. Omit `--no-solvent` only when you intend solvent-box minimization.

## Defaults and limits

- Whole-system minimization is attempted by default, including compatible
  glycans and ligands. Every retained component needs suitable atoms, bonds,
  and force-field parameters.
- Temporary solvent is not retained in the output. `--keep-water` preserves
  crystallographic waters that were already in the input.
- Existing hydrogens are kept when possible. `--rebuild-h` requests rebuilding;
  variant handling may also require a strip-and-rebuild path.
- Protein protonation-variant identities are retained; a separate `protonate`
  stage is not required simply to restore their names.
- Complete chains with identical atom identities and coordinates produce a
  warning. Check for accidentally merged coordinate frames before proceeding.

## Three-Tier Restraint System

| Atom category | Default strength (kcal/mol/Å²) | Purpose |
|---|---|---|
| Original protein heavy atoms | 100 | Keep the observed structure near its starting geometry |
| Rebuilt backbone, including CB | 5 | Restrain repaired regions weakly |
| Heterogen heavy atoms | 5 | Limit ligand/glycan drift |
| New sidechains and hydrogens | 0 | Allow local relaxation |

Two phases use the configured strengths and then tenfold-weaker restraints.
`--max-iter` limits each phase. Keep matching `.dat` metadata so rebuilt atoms
receive the intended tier. Without it, protein heavy atoms are treated as original.

## Usage

### Choose a force field

```bash
dvbfixer minimize prepared.pdb --dat prepared.dat --ff amber \
  --no-solvent -o minimized.pdb
```

The default is `--ff auto`; read its reported choice. GLYCAM-named input selects
compatible AMBER/GLYCAM templates. See [force fields](../force-fields.md) before
supplying custom aliases, XML files, or `--extra-ff` parameters.

### Require ligand parameterization

```bash
dvbfixer minimize prepared.pdb --dat prepared.dat --parametrize-ligands \
  --no-solvent -o minimized.pdb
```

Eligible unknown organic ligands trigger a non-strict automatic GAFF2/AM1-BCC
attempt even without this flag. Explicit `--parametrize-ligands` makes failure
fatal instead of falling back. Requires AmberTools and `openmmforcefields`;
templates are cached under `~/.cache/dvbfixer/lig_params/`.

Complex cofactors such as FAD and HEM are rejected before generic GAFF, even in
automatic non-strict mode. Supply validated state-specific parameters. Parameter
XML does not reconstruct missing atoms. Cross-residue covalent ligands also need
compatible linkage templates.

### Exclude heterogens from optimization

```bash
dvbfixer minimize prepared.pdb --strip-heterogens --no-solvent -o protein_relaxed.pdb
```

Despite its name, `--strip-heterogens` restores heterogens to the output afterward.
It does not jointly relax a noncovalent ligand and its pocket; a restored ligand
can clash with the moved protein. Other automatic setup failures may also reach
this legacy protein-only fallback. Read warnings and inspect the output rather
than assuming a whole-system calculation ran.

### Optional geometry refinement

```bash
dvbfixer minimize prepared.pdb --dat prepared.dat --no-solvent \
  --obminimize-refine --refine-heterogens-only -o refined.pdb
```

OpenBabel UFF refinement is an optional post-pass; xtb GFN-FF is available via
`--xtb-refine`. Neither can bypass an earlier OpenMM setup rejection or guarantee
support for arbitrary chemistry. These are external classical-force-field
geometry optimizations, **not quantum calculations or a dedicated protein-aware
ligand refinement protocol**. In particular, DVBFixer invokes xTB with
`--opt --gfnff`, not GFN1-xTB or GFN2-xTB.

With `--refine-heterogens-only`, most of the protein is **absent from the
calculation**, not merely frozen. The subsystem contains all heterogens and
full covalently linked anchor residues, with selected anchor atoms fixed.
A noncovalent ligand therefore has no surrounding protein pocket in this
post-pass. Its internal strain may decrease while its contacts with the omitted
protein worsen. The preceding OpenMM pocket-aware result does not constrain
this subsequent optimization. Whole-system mode includes the protein but uses
the external force field, not the preceding OpenMM parameters or restraint tiers.
See the [technical scope and limitations](../force-fields.md#technical-scope-of-the-external-post-passes).

UFF is the default OpenBabel force field. MMFF94 variants can mis-type N-linked
sugar geometry. xTB has an upstream fixed-atom drift fix in
[PR #1429](https://github.com/grimme-lab/xtb/pull/1429), merged after the latest
published v6.7.1 release as checked on 2026-10-27. Availability in a particular
package build and correctness on DVBFixer inputs still require verification;
see the [xTB review checklist](../plans/environment-dependency-audit.md#b12-xtb-на-карандаше-проверить-до-решения-об-установочном-профиле).
Inspect linkages rather than assuming a universal-force-field pass improves them.
For systems above 5,000 atoms, whole-system post-refinement can be switched to
heterogens-only with an informational message. Read the run log.

## Options

Frequently adjusted options are `--ff`, `--dat`, `--restraint-k`, `--weak-k`,
`--max-iter`, `--platform`, and `--atom-naming`. Use `--atom-naming standard`
for native output labels or the default `gromacs` for its consumer aliases.
The [generated reference](../reference/minimize.md) lists all defaults and choices.

## Recommended Workflow

Follow the [manual recipe](../pipelines.md#manual-step-by-step) for connected
filenames and sidecars. Do not reuse repair metadata after a viewer silently
changes chain or residue identities. Inspect diagnostics, modeled regions,
ligands, and covalent links after relaxation.

## How it works

Minimization applies two restraint phases and chirality protection. Emergency
chirality repairs are recorded in `REMARK 999 DVBFIXER CHIRALITY_REPAIR`; inspect
their local hydrogen angles. Unresolved D geometry aborts output. Seeded repair
does not imply bit-identical coordinates across numerical platforms.
See [architecture](../../ARCHITECTURE.md) for implementation contracts and
[historical troubleshooting](../known-issues-history.md) for resolved issues.

## Batch mode

```bash
dvbfixer minimize --input-dir prepared --output-dir minimized --no-solvent
```

Options are shared across inputs. See [batch mode](../batch-mode.md).

## See also

- [prepare](prepare.md) — produce the matching structure and sidecar
- [protonate](protonate.md) — explicitly reassign states at another pH
- [diagnose](diagnose.md) — inspect structural findings
- [Troubleshooting](../known-issues.md) — setup failures and fallback limitations
