# New-input geometry qualification: the earlier 12/12 result does not generalize

Status: **bounded prospective input/kernel qualification, not production-ready
coordinate repair**. Continues [source-geometry controls](start-geometry.md).
All initial refusals, generated failures and unfavorable accuracy remain retained.

## Saved-force disaggregation of the earlier experiment

Before interpreting a new cohort, a **retrospective saved-only** audit compared
the existing serialized forces and starting coordinates; it ran no optimizer or
native-reference scoring. [Force evidence](start-geometry-mediators-evidence.json):

- All **eight THR distal-only contrasts** have an exactly identical serialized
  **entire force list** to their raw counterparts. Starting heavy/H positions differ.
- All **eight THR full-rebuild contrasts** change only the generated-CB positional
  anchor and target-Cα guard reference volume. Restraint weights and other
  serialized force terms do not change.
- The **eight VAL no-op contrasts** retain identical controls.

Thus the distal-only effect is not explained by a changed Cα-wall reference or CB
restraint anchor in those cases. H and heavy starting coordinates still change
together; this does not identify their separate contributions, prove equality of
every possible runtime property or explain the original unforced THR2 crossing.
The analysis is explicitly retrospective, not a new pre-registered native trial.

## Pre-download enrollment and source admission

Before downloading any coordinates, `prospective_geometry_cohort.py` froze
**1FNA**, **1SHG** and **1VII**, the previous geometry kernels and all three arms.
The rule remains first coordinate model, first ATOM chain, first two THR and first
two VAL identities, whole-sidechain masks. Source H, altloc, occupancy and touching
LINK/SSBOND/CONECT endpoints stay intact. No source replacement or cleanup is
allowed to gain admission.

Transport/readback code was completed around the enrolled kernels; the independent
auditor was finished after download but hash-bound in the geometry ledger before
native reporting. This is prospective input/kernel policy, not a fully
preimplemented, randomized or blinded benchmark.

All three official RCSB downloads returned coordinate PDBs with retained original
bytes, URL/time and SHA-256 provenance. The source rule produces:

| Source | Masks | Admission outcome |
|---|---:|---|
| 1FNA, X-ray fibronectin III | 4 | All refused by the unchanged heavy-only source scope; 150 deposited H records retained |
| 1SHG, X-ray spectrin SH3 | 4 | All admitted and prepared |
| 1VII, solution-NMR villin headpiece | 2 | Both refused by the unchanged heavy-only scope; 301 deposited H records retained; one THR/one VAL count shortfall |

These are **valid source structures outside this private input capability**, not
invalid PDB depositions. Neither H-containing source was silently dehydrogenated.
Consequently only **one protein, four masks** contributes executed new science;
six correlated source refusals remain visible in the ten-mask denominator.

One legacy prepare runs per admitted mask, with Reference explicitly selected for
PDBFixer and H placement, H seed 1, unchanged PROPKA/Reduce/pH. All four preparations
complete. Refinement reuses the previously frozen `raw`, `beta_only` and
`full_generated` recipes, both H seeds, forced H rebuild and Reference inside H
placement and native minimization. No new angle, force-strength, tolerance, seed
grid or profile choice is introduced.

## Geometry is committed before native accuracy

[Geometry ledger](prospective-geometry-evidence.json) was generated and committed
before parsing held-out reference coordinates. All **24/24 refinement cells**
complete, **48 native phases**, zero execution blockers. The target-geometry/
global-stereo gate is unchanged; initial failures are not used to filter execution.

| Arm | Executed cells | Common research gate passes | THR gate passes |
|---|---:|---:|---:|
| Raw | 8 | **2** | 0/4 |
| Distal-only | 8 | **3** | 1/4 |
| Full generated sidechain | 8 | **5** | 3/4 |

Eight cells per arm are four masks × two H seeds, **not eight independent
proteins**. The previous six-mask 12/12 result is not a guarantee on new inputs.
No refused final structure is promoted into accepted geometry merely because its
target Cβ is R.

### Distinct reasons for rejection

- **VAL23, both seeds/all arms:** the existing gate reports
  `val_branch_name_orientation`. VAL's methyl branches are chemically equivalent;
  this is a stereospecific atom-name convention check, **not an inverted Cβ
  stereocentre**. The frozen gate is retained, not relaxed after seeing results.
- **THR24 in raw/distal-only, both seeds; THR32 raw, both seeds; THR32 distal-only,
  seed 2:** target Cα becomes near-planar and has angle findings. RDKit does not
  assign target Cα CIP in these flat outcomes; target CB remains R.
- **THR24 full-rebuild, seed 1:** the target itself passes, but the independently
  audited *observed context* residue **GLN16 Cα** becomes near-planar. Its retained
  signed volume changes from **2.528406 Å³** before phase 1, to **0.529063 Å³** after
  phase 1 and **0.482378 Å³** after phase 2. The volume remains positive: this is
  **not a D-Cα finding**. The global research gate still refuses the output.

The last result is especially important: fixing generated target geometry cannot
stand in for a final whole-context stereo/shape audit. A finite relative-volume
Cα wall is not an absolute tetrahedral-quality guarantee. This is one retained
case, not proof of its complete force-level cause or a whole-trajectory claim.

## Native accuracy: no hidden removal of geometry failures

[Accuracy evidence](prospective-geometry-native-evidence.json) reports **all 24
completed outputs**, including geometry failures, only after the geometry ledger.
All three unordered arm pairs, all/THR/VAL strata and both all-completed and
common-passing comparisons are retained: **18 descriptive matched rows**. No
native-error or favorable-pair selection chooses a method.

On the **same four THR cells**, original-frame exact-name generated-heavy RMSD
medians are:

| Arm | Median RMSD, Å | THR common gate passes |
|---|---:|---:|
| Raw | 0.509956 | 0/4 |
| Distal-only | 0.519561 | 1/4 |
| Full generated sidechain | 0.927699 | 3/4 |

Full reconstruction improves the geometry pass count here but does **not** improve
the held-out coordinate median. Conversely, the lower raw RMSD does not make its
flat target Cα outcomes acceptable. The common-passing raw/full THR intersection
is empty, explicitly reported rather than assigned a favorable accuracy claim.
Across all eight completed cells, the corresponding exact-name medians are
1.029991/1.043625/1.415851 Å; mixed VAL/THR totals are not a population ranking.

All four admitted native controls pass the unchanged common gate. VAL symmetry
reporting is included separately; no native-driven atom renaming, candidate
selection or observation restoration occurs. Native error is measured in the
original frame without alignment; original-observation motion is retained
separately and does not get reclassified as generated-atom accuracy.

## Qualification outcome and scope boundary

This cohort supports **neither a production repair route nor a general 12/12
guarantee**. It exposes additional near-planar Cα/context failures, a separate VAL
naming limitation and worsening THR coordinate medians. The original spontaneous
THR2 Cβ inversion still is not replayed or uniquely explained. No failed source
or completed output is replaced, silently cleaned or omitted from its denominator.

At the time of this cohort, the production `diagnose` chemistry pass checked Cα handedness; its summary
uses the category `chirality` to mean “D-isomer error”. Adding THR-Cβ or near-planar
shape findings under that category without separating their meaning would be
misleading. A useful next product step would be **detection-only diagnostics with
separate labels**, not automatic coordinate correction or accepting this research
rebuilder as safe. That step changes production behavior and requires explicit
user approval under the standing research-only scope.

That warning-only diagnostic step was subsequently authorized, implemented on a
separate branch, merged into main, and explicitly transferred into this research
branch. Coordinate repair remains unqualified. The subsequent
[saved-force readback](stereo-force-projection.md) explains the limitation of the
relative Cα wall and reports local force competition without claiming original
THR2 causality or changing these frozen geometry/accuracy outcomes.

Private artifacts remain under `tmp/pdbfixer-research/prospective-geometry-`:
pre-download enrollment, originals/provenance, masks, prepared outputs, refinement
plan/run, original geometry ledger/native report and independent saved-only
readbacks. New sources are not registered as tracked production fixtures here.
The separate saved-system analysis uses `start-geometry-mediators-01.json`.

No production/default/dependency/library/public capability, source admission
rules or acceptance thresholds changed. No commit/push. This work uses no
PyRosetta tables; the independent licensing/approval gate for that route remains.

Verification: **825 focused tests passed in 27.85 s**. Targeted Ruff and the
17-file completion-package mypy check pass; all three publications regenerate
byte-identically from saved evidence. This is not full repository/platform acceptance.
