# Troubleshooting and known limitations

[User guide](README.md) · [Installation](installation.md) · [Recipes](pipelines.md)

Start with the symptom, not the implementation history. Keep the full run log
when reporting a problem; warnings often explain why a fallback was selected.

## Find the right check

| Symptom | First action |
|---|---|
| Missing tool, import failure, license error | Run `dvbfixer doctor` in the active environment |
| Structure has clashes, broken geometry, or missing atoms | Run `dvbfixer diagnose input.pdb` |
| GUI cannot find a newly installed executable | Restart its server in the correct environment |
| Model cannot place sequences | Check full sequence and case-sensitive FASTA chain IDs |
| Force-field template mismatch | Check component completeness, naming, and parameters |
| Output naming is unexpected | Check profile and distinguish naming from protonation |

## Save a useful log

Every command accepts `--log-file PATH`. For example:

```bash
dvbfixer prepare input.pdb -o prepared.pdb --log-file prepare.log
```

The log captures Python and inherited child-process output. It appends a versioned
UTC run header and includes a deduplicated warning/error summary when findings
exist. A clean run has no empty warning banner. See [batch mode and logging](batch-mode.md).

## Native Apple installation profiles

Micromamba 2.9.0 and mamba 2.10.0 passed native arm64 Standard/Desktop acceptance.
Micromamba 2.0.2 stalled on the same host; use a recent manager rather than
waiting indefinitely. The helper isolates libmamba profile channels from user
rc additions and clears inherited conda activation state only in child processes,
so a separately selected tools manager cannot accidentally run its own Python
instead of the target prefix's interpreter.

Conda 26.9.4 exposed upstream gettext/libintl library and tk/xorg-libx11 header
path collisions even with a fresh package cache. Strict conflict prevention
correctly rejects that transaction. Do not force `--clobber` to claim a clean
install; use a verified micromamba/mamba route while this native conda gate is
open. An older cache also contained a damaged blosc package: retain evidence and
choose a new cache rather than deleting existing environments or disabling
safety checks.

The declared optional `probe=2.18` package has no osx-arm64 build in bioconda.
Omit `--addon probe` for the native managed route; the existing Python diagnostic
fallback remains available. This does not affect Standard/Desktop or the
independent `msa-alternatives` add-on. A source-built native Probe route has not
been established by this acceptance increment. Do not silently substitute
another Probe version or a Rosetta environment. Actual commands and evidence
are in [profile acceptance](installation-profile-acceptance.md).

## Diagnostics report errors after a completed run

### Diffusion has no independently validated candidate

This means the rejected candidate is not safe for downstream preparation; it is
not permission to bypass geometry or chirality gates. New Apple requests select
`protpardelle-1c-cc94-mps-v1`; its native MPS acceptance remains pending.
Historical cc89 launchers are for explicit frozen-request replay and do not
qualify cc94. See the [cc94 deployment guide](../deploy/protpardelle-1c/cc94.md).

Scientific refinement failures preserve the raw candidate with the specific
exception in its provenance. Runner/validation failures retain a private workspace
under `--diffusion-work-parent`, announced in the log, for operator inspection.
It may contain the full input structure and should not be shared publicly without
review. Delete retained diagnostic directories when no longer needed. GUI/container
hosts can impose their own work-directory lifecycle and quotas.
The [supplied 5J13 regression](research/multi-gap-diffusion-5j13.md) documents
the corrected multi-gap run and its remaining, correctly rejected geometry.

### Ordinary diagnostic findings

`diagnose` exit code `1` means analysis found ERROR-level structural findings.
It is not by itself an execution failure. Inspect the affected residues and
choose the appropriate repair tool; do not just suppress the report.

ZBS warns about postflight findings by default. Use `--strict-postflight` when
those findings must fail a pipeline. Read [diagnose](commands/diagnose.md) for
severity meanings and [zbs](commands/zbs.md) for postflight options.

## Missing MODELLER, Reduce, or sequence aligners

- Configure the [MODELLER license](installation.md#modeller-installation-and-license)
  for gap rebuilding and comparative modeling.
- Install [AmberTools and Reduce](installation.md#ambertools-and-molprobity-reduce).
  Missing Reduce is fatal for explicit `tleap-reduce`; legacy preparation can
  warn and use its documented fallback, which is not evidence that Reduce ran.
- Check [MSA engines](installation.md#multiple-sequence-alignment-executables)
  with `dvbfixer msa --list-engines`. MUSCLE must be v5.

## Antibody numbering tool failed

An explicitly selected antibody scheme requires the ANARCI Python API, `hmmscan`
and its pressed HMM database. Missing data or a tool failure is an error before
output, not "no domain found" and not permission to silently substitute SEQRES.
The same distinction applies to Homology's antibody mode.

## Trajectory/topology or TPR reader error

`cluster` checks the topology/trajectory atom counts and first frame before
clustering. Use inputs with matching atom ordering. If the installed MDAnalysis
cannot read a TPR writer version, supply a supported matching TPR or an explicitly
chosen matching GRO/PDB topology. This is not an automatic topology substitution.
Requested Plotly output is checked before analysis; its documented warning and
numerical-output fallback are retained when the optional plotting module is absent.

## Missing sequence or incorrect chain mapping

Modeling needs the intended complete sequence, from SEQRES or `--fasta`.
FASTA headers must encode chain IDs: `>chain_H`, `>1abc_H`, or `>H`.
Chains `H` and `h` are different; sequences are not assigned by file order.

Use `--no-terminal` to omit missing N/C tails, not to erase internal gaps.
Inspect repeated-sequence alignments and modeled regions. See [model](commands/model.md).

## No force-field template for a ligand or cofactor

Naming, atom completeness, and parameterization are separate requirements.
Supplying an XML parameter file does not reconstruct missing cofactor atoms.

- Use [component-info](commands/component-info.md) to inspect CCD chemistry and
  heavy-atom completeness.
- Use authoritative `prepare --smiles` only for matching isolated single-residue
  ligands; it does not support externally covalent components.
- Minimize automatically attempts GAFF2 for eligible unknown organic ligands.
  Explicit `--parametrize-ligands` makes parameterization failure fatal.
- Protected cofactors such as FAD/HEM require validated state-specific parameters;
  they are not admitted to generic GAFF merely because automatic mode is non-strict.
- `minimize --strip-heterogens` excludes heterogens from optimization and restores
  them afterward. It does not jointly relax a noncovalent ligand and its pocket.

Inspect warnings and retained components in the output. See [minimize](commands/minimize.md)
and [force fields](force-fields.md).

## Glycan naming, connectivity, or drift

Use the [GLYCAM recipe](pipelines.md#glycoprotein-preparation-with-glycam) for
compatible AMBER/GLYCAM preparation and ACPYPE export. PDB names such as NAG do
not by themselves select the corresponding GLYCAM templates.

Preparation uses NLN/OLS/OLT only when the attached sugar is GLYCAM-named.
Automatic CONECT inference supplements input connectivity, but inspect links
and sugar identities rather than assuming all deposited annotations are complete.
Use [conect](commands/conect.md) for an explicit preprocessing report.

Glycan heavy atoms have weak positional restraints during minimize, not zero
restraints. A grossly displaced or clashing tree still requires investigation.
Optional universal-force-field post-passes do not replace valid initial setup.

## Unexpected residue or hydrogen names

- `rename` folds variants to parent names without rebuilding H.
- `atom-names` translates target-specific labels and preserves known variants by
  default. Its opt-in inference reads represented H patterns, not predicted pKa.
- `protonate` chooses supported states and rebuilds hydrogens.

Choose the intended [naming tool](pipelines.md#choose-a-naming-tool). Use native
`standard` or consumer-specific `gromacs` naming deliberately; they are not
interchangeable. A separate protonation run is not needed just to restore labels.

AMBER14 lacks terminal protonated ASP/GLU templates. A terminal ASH/GLH can be
reverted to ASP/GLU with a warning on affected paths. Do not interpret this as
successful preservation of the requested terminal charge state.

## Sidecar no longer matches the structure

Keep model/prepare `.dat` files beside their matching PDB. Saving through a viewer
that reassigns chain IDs or residue numbers can invalidate those identities.
Return to matching artifacts or regenerate preparation metadata rather than
reusing an unrelated sidecar.

## Large CIF structures exceed PDB limits

CIF input is normalized to internal fixed-column PDB. At most 62 chain IDs are
representable; atom/residue identifiers and coordinates must also fit. CIF
normalization refuses unrepresentable data rather than truncating it. Do not
assume downstream PDB CONECT writers support six-digit atom serials.
Split the system or use a format-capable external workflow when necessary.
See [CIF input](cif-input.md).

## GROMACS water, ions, or mixed scaling

AMBER RTP export matches ion parameters to `--water` by default; do not mix
unrelated ion and water sets. CHARMM only accepts its supported water choices.
ACPYPE uses TIP3P and ignores `--water`/`--ion-set`.

AMBER and GLYCAM have different 1–4 scaling. Use `top --acpype` for that combined
route rather than assuming one global GROMACS fudge factor covers both.
If an ACPYPE `interchain_ss.itp` include contains intermolecular interactions,
keep it below all molecule entries after solvation and ion addition.
See [top](commands/top.md).

## Chirality warnings and reproducibility

Inspect residues marked `REMARK 999 DVBFIXER CHIRALITY_REPAIR`, especially local
hydrogen angles. Minimize guards initially-L centers and rejects unresolved
D geometry; rerunning until a lucky result appears is not a repair strategy.
Seeded heavy-atom rebuilds do not imply bit-identical whole-system coordinates
across CPU/GPU backends or thread configurations.

## Experimental diffusion limits

### Alternate locations

Diffusion requires one state per residue and no altloc labels, including residues
far from the modeled gap. Use [Altlocs](commands/altlocs.md) in the GUI or CLI to
review/select states and write a separate single-state PDB. Input checks reject
unresolved altlocs before sampler launch; they are not silently resolved by atom
encounter order. The final validator independently checks the candidate too.

`external-covalent-link-near-generated-region` refers to unsupported
linkage chemistry involving non-fixed atoms near a generated region. CONECT
bonds within a fixed anchor residue and links between two fixed atoms (including
disulfides) are allowed even at gap boundaries. The source coordinates and link
records must be preserved, and independent validation checks them. If a refusal
remains, the error lists exact `chain:residue:atom` endpoints. Inspect the masks
and actual chemistry; do not delete a genuine link to bypass a mismatch.

Gap length outside the selected profile's validated range is advisory: the tool
warns and attempts reconstruction. This includes internal gaps shorter than three
residues and long terminal tails. Use `--no-terminal` to crop tails instead.
This does not waive independent output validation or backend-specific topology,
chemistry, context, or resource restrictions.

Diffusion requires a separately installed compatible runner and checkpoint.
Availability, gap sizes, terminal support, ligands, and batch/GUI integration
depend on the selected profile and workflow. There is no silent MODELLER fallback.
Start with the current [model scope](commands/model.md#experimental-diffusion-scope)
and [runner installation](installation.md#experimental-diffusion-runner-environments),
not historical acceptance notes. A passing `doctor` handshake does not establish
broad scientific validation.

## Environment and developer notes

- Missing `libmkl_gf_lp64.so.3`: follow the [OpenBLAS repair](installation.md#numpyscipy-blas-backend).
- PROPKA on Python 3.14: retain the supported Python `<3.14` constraint.
- NumPy typing errors: development uses NumPy `<2.5` with a Python 3.11 baseline.
- macOS Docker ncurses symlink failures: put the environment/cache on the
  container's native filesystem, not a case-insensitive host bind mount.
- OpenMM CUDA/PTX errors: verify host-driver/toolkit compatibility before
  changing scientific settings.

Detailed environment diagnoses and resolved bugs are retained in the
[historical troubleshooting archive](known-issues-history.md). That archive is
implementation history, not the recommended starting path for users.
