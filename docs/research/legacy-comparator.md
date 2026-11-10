# Full legacy comparator: real observation movement and two THR β-stereo failures

Status: **bounded retrospective comparator on six previously admitted masks**,
not a production change, general method ranking or independently qualified cohort.
Continues the [prospective six-arm comparison](prospective-chi-guard.md), whose
bare seeded PDBFixer arm was not the full production pipeline.

## Frozen execution

Before execution, `research_legacy_comparator.py` froze all **six shared-admitted
new-cohort masks**, including the terminal THR unsupported by library queries.
No native-error, stock-success or geometry-success filtering. The four prior
Top7 projection refusals remain outside shared admission, not silently repaired.

Each case ran once, in a fresh owned directory, through the unchanged CLI:

```bash
dvbfixer prepare input.pdb --backend legacy --atom-naming standard \
  -o prepared.pdb --dat prepared.dat
dvbfixer minimize prepared.pdb --no-solvent --atom-naming standard \
  --dat prepared.dat -o minimized.pdb
```

The runner uses the already established scientific interpreter and source tree,
prepends its bin directory in the child PATH, and preserves every raw output,
sidecar, log and exit status. PROPKA/Reduce, default pH/force field/platform,
original/weak restraint strengths and iteration budgets remain unchanged.
`--no-solvent` is the explicit development convention, not a new optimization
profile; standard atom naming is an output/audit convention. This is not ZBS or
Modeller gap filling. No external installation or production file was modified.

All **6/6 prepare and minimize workflows completed**. Completion does not imply
the common research geometry gate passed. No repeatability claim follows from
one automatic-platform execution per case.

## Geometry is measured in the actual output frame

The auditor never puts additions back onto original observations. It creates a
clearly labelled, private **heavy-only actual-frame diagnostic view**, remasking
only the declared target atoms for the existing common validator. Every retained
heavy coordinate stays exactly at its real output position. Shared variant-name
normalization is applied only to a private copy; raw H-containing outputs remain
unchanged. H geometry is **not** covered by this common heavy-atom audit.

The validator's `observed_preserved` field refers only to this actual-frame
diagnostic view. **It does not mean legacy preserved original observations.**
Original-input heavy identities and displacement are measured separately, without
alignment, including explicit missing-identity accounting. All original heavy
identities were retained here. Preparation preserved their coordinates; minimize
moved all of them at PDB precision in each case.

## Results

[Evidence](legacy-comparator-evidence.json); generated-target exact-name RMSD in
the original reference frame, Å. No superposition or post-hoc coordinate repair:

| Mask | Prepare RMSD | Minimize RMSD | Final common geometry | Maximum observed movement |
|---|---:|---:|---|---:|
| 1ENH THR6 | 1.996717 | 1.946369 | Reject: CB S | 1.306426 |
| 1ENH THR27 | 3.110846 | 1.910700 | Pass | 0.969650 |
| 1CRN THR1 | 1.834284 | 0.817626 | Pass | 0.710230 |
| 1CRN THR2 | 2.652091 | 1.509790 | Reject: CB S | 0.736986 |
| 1CRN VAL8 | 0.262753 | 0.270669 | Pass | 0.741395 |
| 1CRN VAL15 | 0.072146 | 0.186862 | Pass | 0.749131 |

- Prepare common gate: **2/6**. Three THR Cα centres have small positive signed
  volumes below the unchanged research near-planarity bound; these are **not
  negative-volume D-Cα**. Other angle/clash findings remain in the retained audit.
- Minimize common gate: **4/6**. All target Cα centres are S, with positive volumes
  about 2.45–2.58 Å³. The two final failures concern **THR CB**, not D-Cα.
- `1ENH THR6` has CB S after both prepare and minimize. `1CRN THR2` has CB R after
  prepare and **CB S after minimize**. Independent target RDKit CIP readback and
  named β-volume checks agree. This identifies a sidechain stereochemistry issue
  separate from the Cα guard; it does not establish the complete causal mechanism.
- Original-heavy RMS displacement after minimize ranges **0.120–0.216 Å**;
  maximum individual displacement **0.710–1.306 Å**. No observation was restored.

Across all six masks, prepare target RMSD median is **1.915500 Å** and minimize
median **1.163708 Å**. On these same six masks the previously committed synthetic
baseline median is **0.191971 Å**, with six selected common-pass outputs and fixed
observations. These deliberately selected, correlated cases from two proteins do
not establish broad superiority, population accuracy or a production replacement.
The earlier **stock-supported matched-two accuracy counterexample remains valid**;
the full legacy workflow is a different comparator, with different movement and
outcomes. No arm was reranked or tuned using these native errors.

## Reporting correction retained explicitly

The initial saved legacy report reused a fixed-observation reporting helper. Its
target RMSDs were valid, but its auxiliary χ error combined original N/CA anchors
with relaxed additions. A separate hash-bound **actual-frame χ reporting
amendment** corrects only that field from the complete real output. Original v1
report and auditor remain retained; **all RMSD, motion, geometry, coverage and
choices are unchanged**. The published v2 report binds the amendment and superseded
report hashes. This was a reporting correction after native metrics were seen,
not a prospective scientific change or another sampling run.

## Retained evidence and next question

All private artifacts are under `tmp/pdbfixer-research/`:
`legacy-comparator-plan-01.json`, `legacy-comparator-run-01/`,
`legacy-comparator-audit-01/selection-before-native.json`,
`legacy-comparator-report-01.json`, and `legacy-chi-amendment-plan-01.json`.
`research_legacy_comparator.py` binds production source hashes and inputs;
`audit_legacy_comparator.py` independently reads saved outputs and records geometry
before native reporting. `amend_legacy_chi_reporting.py` performs the explicit
saved-only χ correction. Inputs and every output/log digest remain verified.

The next meaningful mechanism question is why **THR2's β centre changes during
legacy minimization**, and how to protect generated sidechain stereocentres
without moving observations or silently swapping chemically distinct THR branches.
That requires a separately declared diagnostic experiment, not an immediate
production patch, adaptive optimizer grid or reinterpretation as a Cα failure.

No production/default/dependency/public completion capability changed; no
commit/push. The installed PyRosetta route remains research-only with the
[original-table/provenance/license/explicit approval gate](../plans/pdbfixer-reconstruction-handoff.md#current-authorization-installed-pyrosetta-for-research-only-2026-11-03)
unchanged. Full legacy execution itself does not depend on that library.

Verification: **800 focused tests passed in 26.25 s**. Saved-only audit ledger,
initial report and corrected publication regenerate byte-identically. This is not
the full repository/platform acceptance suite or a production qualification.

A subsequent [logging-only native phase trace](legacy-stereo-trace.md) does not
reproduce the original endpoints: THR2 stays CB R in that replay, while THR1 develops
a near-planar Cα. The original failures remain retained; no phase localization or
causal claim follows from a nonexact replay.
