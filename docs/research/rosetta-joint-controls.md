# Exact two-target scoring controls: exhaustive search without a packing claim

Status: **bounded joint-mechanism controls**, not new independent-family accuracy
qualification, a public completion backend or an interaction-graph approximation.
Follows the [predeclared single-target follow-up](rosetta-scoring-followup.md).

**Main result: this cohort does not demonstrate a benefit from coupled packing.**
Exact search and multi-target preservation work, but almost all candidate matrices
have only one varying axis. The sole genuine 3×3 score matrix is numerically
additive on the tested grid. Improved selections must not be called a joint
interaction-optimization success.

## Enrollment and bounded search

Before new joint reconstruction/scoring, `scripts/research_rosetta_joint.py enroll`
froze **12 cases** using already studied deposited 1CRN and 1AKI sources. Select
the first two VAL and first two THR identities in chain-A file order, group their
file-order list into two consecutive pairs, and apply three mask classes:

| Group | Residue pair | Mask classes |
|---|---|---|
| 1CRN group 0 | THR1 / THR2 | Distal/distal; whole/whole; first anchored/second distal |
| 1CRN group 1 | VAL8 / VAL15 | Same three classes |
| 1AKI group 0 | VAL2 / VAL29 | Same three classes |
| 1AKI group 1 | THR40 / THR43 | Same three classes |

"Whole" removes CB and both distal atoms. "Anchored" removes only the first named
distal atom (OG1 for THR, CG1 for VAL), retaining CB and its partner. Original
observations/link records survive unchanged; no H stripping or repair. Every case
is retained, including the known difficult VAL8 anchored constituent.

The systems and constituent single-target outcomes were already known. New joint
combinations were enrolled before execution, but this is **not** an independent
prospective accuracy cohort, randomized sample or new protein-family validation.
The scope deliberately excludes new chemistry and production consumers.

The existing `plan_mixed_cohort` / `complete_mixed_cohort` kernels are unchanged:

1. Generate and audit each target against the shared masked source.
2. Retain all attempts, including source-only rejections.
3. Enumerate the complete Cartesian product of **individually feasible** pools.
4. Audit each joint combination with the other target's generated atoms present.
5. Apply canonical rounded materialization and global Cα/THR-Cβ guards.
6. Require independent target RDKit CIP before energetic eligibility.

At two targets and three proposals each, this is at most **nine combinations per
case**. Exactness refers to this restricted pool, not all conformations or a global
physical optimum. The original four-target/81-combination kernel limits are intact.

## Fixed-atom scoring and independent readback

The frozen single-target `score_cell` is called on the **same full candidate twice**,
once for each target graph. Both imports must preserve all physical heavy input
identities/coordinates exactly, match external links and independently match their
own target's canonical internal bonds. Both views must have identical full heavy,
H, virtual, residue/type, external-link and allowed-extra snapshots and equal finite
scores. Neither a target omission, score average nor near-equal substituted pose is
accepted. This redundant validation is a mechanism control, not an optimized runner.

Within each case, residue/type/topology inventories and fixed-context H, terminal
OXT and virtual coordinates must agree across candidates. Only target H may vary
with the proposed heavy coordinates. `ref2015`, initialization flags, build
`1a56185c2592611dec4c9c75ddc9468cd2227c1f`, single threading and seed 7 are unchanged.
No packer, minimizer, reflection or post-hoc observation restoration is called.
Context internal bond orders/default protonation are not independently qualified;
the H protocol is not a replacement for legacy PROPKA/Reduce. Runtime build and
weights are recorded, not a new cross-machine/full-database qualification.

Select minimum full-pose score only within the **complete old raw/rounded joint
feasible set**; an unscored feasible combination blocks selection. Commit before
native imports. Canonical generated-only `DatRecord` sidecars are made during
materialization; private selected PDB/`.dat` copies must match the winning cell.
No Rosetta-created H/OXT/virtual atom enters added-atom ownership.

`scripts/audit_rosetta_joint.py` replays source planning, individual audits/winners,
the exact Cartesian enumeration, joint audits, rounded output and canonical `.dat`
from **saved proposals**, without calling a coordinate generator. It independently
checks both graph views, exact snapshots, score-term sums/weights, fixed context,
global and target stereo, source records and the selected private artifact bytes.

[Portable independent evidence](rosetta-joint-controls-evidence.json) contains
per-target and pooled held-out errors and the enrolled coupling diagnostic.
Pooled RMSD is over all generated atoms across both targets, not a mean of two
residue RMSDs. VAL swap reporting is allowed only for both-generated methyls;
fixed partners and THR branches never swap. Truth is reporting only.

## Executed counts and limitations

| Measure | Count |
|---|---:|
| Enrolled / source-admitted joint cases | 12 / 12 |
| Complete source-feasible Cartesian combinations | 34 |
| Joint raw / rounded rejections among these combinations | 0 / 0 |
| Verified generated-only sidecars | 34 |
| Native joint controls, never eligible for selection | 12 |
| Scored full poses / target-view imports | 46 / 92 |
| Cases with selected private structures | 11 |
| Empty source-feasible matrix | 1 |

All 46 full poses pass physical/topology snapshot checks, global Cα/THR-Cβ and
independent target stereo readback. The empty case is **1CRN VAL8/VAL15 anchored+
distal**: VAL8 CG1 has no source-feasible attempt. Native controls do not rescue it,
and no selected artifact is written. Both targets and all reference records remain
retained. The existing chirality/steric gates and historical outputs are unchanged.

No independently selected per-target clearance combination violates the joint
guards in this cohort. Hence it also fails to scientifically challenge the
new–new clash advantage; unit regressions test that mechanism separately.

Compared with the existing **raw-joint worst-clearance selector**, score changes
two choices. Rounded pooled exact-name RMSD:

| Pair/mask | Raw-joint clearance | Exact score choice |
|---|---:|---:|
| 1CRN THR1/THR2, distal/distal | 1.739 Å | 0.183 Å |
| 1AKI VAL2/VAL29, distal/distal | 1.707 Å | 0.114 Å |

Same-case medians across 11 feasible cases: **0.311→0.211 Å**. Do not compare these
pooled two-target numbers with earlier single-target/mutual-template denominators.
This is changed candidate ranking, not increased feasible coverage.

For 1CRN, the score winner already equals the **independent per-target clearance**
combination. For 1AKI, only VAL2 changes (2.414→0.147 Å), repeating the previously
observed single-target effect; VAL29 remains 0.066 Å. Neither is new evidence of
coupled energetic selection or an independent accuracy success.

## Nonadditivity diagnostic: principal negative result

Before execution, the diagnostic was fixed as:

```text
Delta(i,j) = E(i,j) - E(i,j0) - E(i0,j) + E(i0,j0)
```

Reference indices are the first source-feasible indices, never native or best
held-out values. The complete Cartesian score table is required; missing scores
invalidate the diagnostic. This measures conditional **full-score** nonadditivity
in Rosetta units, not a calibrated atom-pair energy or kcal/mol interaction.

- Ten matrices have a single-option axis: their zero delta is **uninformative**,
  not proof of absent physical interactions.
- One matrix is empty and supplies no coupling measurement.
- Only **1AKI VAL2/VAL29 distal/distal** is 3×3. Its observed Cα distance is
  13.128 Å; maximum absolute delta is **1.7053×10⁻¹³**, numerical cancellation scale.
  The tested score table is effectively additive under this fixed protocol/grid.

Thus this cohort validates exhaustive machinery, not the scientific necessity of
an interaction graph. Do not implement or qualify approximate joint packing from
these results. A closer-contact, predeclared cohort with retained zero/one-option
pools is needed; never filter pairs by reconstruction success to manufacture an
apparently coupled benchmark. Richer feasible sampling remains a separate gate.

## Reproduction and scope

Retained local evidence: `tmp/pdbfixer-research/joint-enrollment-01/`,
`joint-plan-01/`, `joint-score-01/` and `joint-audit-01.json` (initial analysis).
`joint-score-01/<case>/selected-candidate.{pdb,dat}` are private research artifacts,
not public diffusion/completion bundles or GUI-registered outputs.

```bash
# Normal scientific environment; prior frozen single-target policy is required.
PYTHONPATH=src:. python scripts/research_rosetta_joint.py enroll --output tmp/pdbfixer-research/new-joint-enrollment
PYTHONPATH=src:. python scripts/research_rosetta_joint.py prepare --enrollment tmp/pdbfixer-research/new-joint-enrollment --output tmp/pdbfixer-research/new-joint-plan
# Existing licensed PyRosetta environment; use exclusive new destinations.
PYTHONPATH=. python scripts/research_rosetta_joint.py score --plan tmp/pdbfixer-research/new-joint-plan --output tmp/pdbfixer-research/new-joint-scores
PYTHONPATH=src:. python scripts/audit_rosetta_joint.py --plan tmp/pdbfixer-research/new-joint-plan --scores tmp/pdbfixer-research/new-joint-scores --output tmp/pdbfixer-research/new-joint-audit.json
```

**611 focused tests passed in 24.97 s**, targeted Ruff passed. Tests prohibit
generation during enrollment/readback, reject incomplete matrices and altered
candidate audits, enforce both target views, and distinguish a degenerate axis
from an informative coupling measurement. This is not a full repository suite,
independent energy-function implementation or broader chemistry qualification.
Canonical saved-candidate readback is byte-identical; package mypy, generated
CLI/GUI checks, agent maps, old evidence and all 49 fixture hashes also pass.
No production changes, new dependency, empirical table export or approximate
interaction graph is introduced.
