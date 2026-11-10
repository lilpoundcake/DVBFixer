# Saved-force readback: relative Cα protection is not a shape guarantee

Status: **retrospective static mechanism analysis**, not coordinate repair,
optimizer replay, native-accuracy qualification or resolution of the original
THR2 Cβ inversion. Continues the [new-input geometry cohort](prospective-geometry.md).
[Complete evidence](stereo-force-projection-evidence.json) retains every enrolled
cell and measured focal centre, including unfavorable states.

## Diagnostic transfer and preserved scientific baseline

Before continuing, the user required the new warning-only diagnostics on this
research branch. Main's `039e1ba` and `9e71372` were cherry-picked as `471844a`
and `8eeeada`, without a full main merge. All 324 pre-existing modified/untracked
research files were checked unchanged immediately after transfer. The active
diagnose/runtime suite passed **164 tests**. The new pure stereocentre measurement
is used here for angle readout; see its [criteria](stereo-diagnostic-criteria.md).

Earlier scientific plans intentionally bind the old production-source hashes.
Those hashes and publications were **not rewritten**. A detached `4a328a6`
historical worktree plus checked copies of the original research files is retained
at `tmp/pdbfixer-research/pre-stereo-baseline-01/`; it shares the original evidence
through a read-only-use symlink. All **30 historical final verifiers** pass from
that original source context. New analysis binds the active post-transfer source
and explicitly verifies prior provenance against the preserved historical root.
This is not a claim that old frozen requests silently run against new source.

The first active focused run correctly exposed one source-inventory mismatch
(838 passed, one failure): its retained diagnostic-consumer file hash was old.
An explicit [current-source amendment](pdbfixer-usage-inventory-stereo-diagnostics.json)
now pins the new hash while requiring unchanged imports, call records and totals
(11 constructors, nine seeded rebuilds, 27 method sites). The original inventory,
offline template probe and frozen scientific plans remain unchanged; no broad
skip, old-result hash substitution or new template qualification was introduced.

## Frozen analysis scope

The outcomes and focal GLN16 failure were already known: this is **retrospective**,
not a prospective/blinded cohort. Before force evaluation, the analysis froze:

- ALL **24 completed cells**, two sides of both saved phases: **96 states**.
- Every target Cα, THR target Cβ only, and the previously identified observed
  GLN16 Cα: **240 centre-state measurements** (192 Cα, 48 THR-Cβ).
- Reference-platform **static** evaluations of the saved systems/coordinates;
  no integration, minimization, constraint projection, H placement, extra seed,
  force-strength grid, coordinate restoration or native-reference scoring.
- Separate force-group measurements on cloned systems, with exact serialized
  restoration proving that only measurement-group metadata changed. Group sums
  reproduce total energies and projected forces; each Cα wall projection also
  matches its independent analytic derivative.

Recorded initial global parameters must match every phase system's defaults;
the bound original minimize pipeline has no `setParameter` calls. System XML
must remain identical across each captured native minimization call. All 24 cells
were evaluated successfully; **zero readback blockers**. The independent saved-only
auditor checks complete enrollment, arithmetic, wall activation and phase switches.

For ordered centre X and neighbors A/B/C, use scalar triple product
`v=((A-X)×(B-X))·(C-X)`. The reported force is
`F·∇v / |∇v|`, in **kJ/mol/nm**, over the four Cartesian atom coordinates.
Positive points toward increasing ordered signed volume; negative points toward
decreasing volume. This is **not** a mass-weighted acceleration, constrained
optimizer direction, integrated work or proof of an unrecorded trajectory.

## Why an extremely flat target can leave its Cα wall inactive

The existing wall uses `q=v/vref`, with `vref` captured from starting L geometry:
`E=0.5*k*max(0.2-q,0)^2`, `k=10000 kJ/mol`. It is a finite relative-volume
penalty, not an absolute bond-angle/tetrahedral-quality boundary.

Examples below are **unrounded phase-2 endpoint** measurements, seed 1:

| Target / arm | Starting vref, Å³ | Wall onset, Å³ | Final v, Å³ | Final angle | Wall active |
|---|---:|---:|---:|---:|---|
| THR24 / raw | 0.025934 | 0.005187 | 0.008107 | +0.152° | No |
| THR24 / full generated | 2.555116 | 0.511023 | 2.519380 | +51.669° | No |
| THR32 / raw | 0.327356 | 0.065471 | 0.065190 | +1.235° | Yes |
| THR32 / full generated | 2.698235 | 0.539647 | 2.444482 | +49.085° | No |

Raw THR24 starts almost planar. Its endpoint is still above 20% of that tiny
starting volume (`q=0.312608`), so the wall correctly evaluates to zero even
though the new diagnostic angle screen warns. THR32 can reach a low-angle endpoint
while sitting close to its similarly small relative floor. This explains a
specific limitation of the wall's stated contract, **not** a faulty force sign,
a missing wall, general repair success or an explanation of THR2 Cβ crossing.

Across the full readback, four of the 192 Cα centre-states have an active wall;
24 have absolute out-of-plane angle ≤10°. These include initial poor geometry
and duplicated phase-boundary coordinates; they are **not independent failure
rates or population statistics**. Every state remains in the evidence.

## GLN16: a restoring wall competes with a flattening angle component

For THR24 / full-generated / seed 1, observed-context GLN16 has:

| Saved state | Volume, Å³ | Angle | q | Wall |
|---|---:|---:|---:|---|
| Phase 1 before | 2.528406 | +51.096° | 1.000000 | Inactive |
| Phase 1 after / phase 2 before | 0.529063 | +10.305° | 0.209248 | Inactive |
| Phase 2 after | 0.482378 | +8.964° | 0.190783 | Active |

At its final saved state, projections are:

| Force component | Projection, kJ/mol/nm |
|---|---:|
| Cα wall | **+2487.696** |
| Harmonic angles | **−2427.663** |
| Nonbonded | −344.256 |
| Positional restraints | +135.876 |
| Harmonic bonds | +171.147 |
| Periodic torsions | −32.473 |
| CMAP | +77.259 |
| Net | +67.586 |

The active wall is restoring, and its measured projection matches the analytic
derivative. Its energy is 0.424738 kJ/mol; a small energy does **not** mean a small
local derivative. The largest opposing group here is the harmonic-angle force,
not the Cα wall. That group includes both heavy-atom and hydrogen angles: this
readback does **not** isolate a particular angle or establish H-only causality.
The positive net projection is not a convergence certificate; constraints and
optimizer stopping/path information are outside this static measurement.

The volume stays positive throughout these recorded boundaries: **not a D-Cα
finding**. The phase endpoint is unrounded; the retained final PDB instead gives
8.910° / 0.479232 Å³. Neither representation may be substituted for the other
without stating the rounding boundary.

## The phase switch is measured, not inferred from different coordinates

For ALL **24 cells**, phase-1-after and phase-2-before coordinate snapshots are
byte-identical. At those identical coordinates, non-positional force-group
energies and focal projections match; positional contributions scale by exactly
one tenth within the audit tolerances.

For the GLN16 case above, its positional restoring projection changes from
**+1350.708** to **+135.071 kJ/mol/nm**, while the harmonic-angle projection remains
**−1503.349** and the Cα wall remains inactive. Net Cartesian projection switches
from **+16.654** to **−1198.983**. The local change is exactly accounted for by
the declared restraint softening at that saved configuration. This is useful
mechanistic evidence, **not** an intervention replay establishing the complete
later trajectory or the sole cause of the final defect.

## Qualification and next boundary

The earlier geometry/RMSD failures and H-source admission refusals remain intact.
No gate, seed, candidate or atom-name convention is relaxed. Nothing here
qualifies automatic reflection, branch swapping, a stronger production wall,
observed-coordinate repair or a new completion backend. The original unforced
1CRN THR2 event remains unresolved and is not part of this 1SHG force readback.

The next bounded mechanism question is which angle terms and H/context changes
mediate the large harmonic-angle projections, with complete retained-state
reporting rather than native-error/profile search. That question is separate
from changing production restraints or claiming a qualified reconstruction route.
No PyRosetta tables are read or redistributed; its separate production
provenance/license/approval gate remains unchanged.

Code: `scripts/audit_stereo_force_projection.py`,
`scripts/audit_stereo_force_readback.py`; private enrollment and raw readback:
`tmp/pdbfixer-research/stereo-force-plan-01.json` and
`tmp/pdbfixer-research/stereo-force-readback-01/`. The preserved historical cwd is
for old verification only; new analysis runs from the active research branch.

Both complete static evaluations reproduce byte-identical reports, and the
independent audit publication regenerates byte-identically. These are repeated
evaluations of the same saved states, not additional biological/seed replicates.

Verification: **840 active focused research tests passed in 29.03 s**, plus the
**164 active diagnose/runtime tests** used for transfer acceptance. All **30
historical final verifiers** pass from the preserved baseline, and the new force
checkpoint verifier passes. Targeted Ruff, 17-file completion-library mypy,
generated CLI/GUI checks, agent maps, whitespace and 1,277 local Markdown targets
pass; seven planned documentation files are complete. The first inventory
migration failure is retained separately rather than relabelled as a pass.
This is not a full-repository/slow-suite/platform or repair qualification.

The subsequent [complete angle-term decomposition](stereo-angle-terms.md) now
separates the large GLN16 harmonic-angle projection: at the known final state,
H-containing terms contribute −3793.179 and heavy-only terms +1365.516
kJ/mol/nm. Initial local heavy coordinates and angle parameters match across all
24 contexts, while generated HA geometry differs. This is stronger local
attribution, not yet a controlled H-start trajectory replay or production repair.
