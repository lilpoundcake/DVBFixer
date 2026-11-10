# Backbone-selected close-contact reconstruction controls

Status: **private mechanism research**. Follows the
[initial exact joint controls](rosetta-joint-controls.md), whose only two-variable
matrix was effectively additive. Production, generators, geometry thresholds,
single-target scorer and original joint runner/auditor are unchanged.

**Result:** joint checks reject three genuine new–new clashes, but no case has
two varying feasible axes. Coupled energetic-packing benefit remains unmeasured.
The independent clearance winner remains jointly feasible in every case.

## Frozen input-only enrollment

`scripts/research_rosetta_contacts.py` enrolled pairs before new reconstruction or
scoring, using the same unchanged first-model chain-A projections of deposited
1CRN and 1AKI. Only observed **Cα coordinates**, canonical VAL/THR identities and
file-order positions enter pair ranking; hidden sidechain coordinates, candidate
diversity, occupancy filtering and reconstruction outcomes do not.

Fixed rule: same-chain Cα distance **≤8 Å**, at least three positions apart in
residue file order, ascending distance with file-order pair ties; at most four
pairs per structure, with shared residues allowed. Apply the original three mask
classes (distal/distal, whole/whole, anchored-first/distal-second) to every pair.
Keep empty/one-option pools and never replace unsuccessful cases or structures.
Full eligible rankings and source/code hashes are frozen in the enrollment.

**1CRN has zero eligible pairs**, retained as a selector outcome, not replaced or
recovered by expanding the threshold. All four 1AKI pairs were selected:

| Rank | Pair | Observed Cα distance |
|---|---|---:|
| 0 | THR89 / VAL92 | 5.169 Å |
| 1 | VAL2 / THR40 | 6.355 Å |
| 2 | THR43 / THR51 | 6.549 Å |
| 3 | THR51 / THR69 | 7.972 Å |

This is closer **backbone** geometry, not a guarantee of contacting sidechains.
Systems and some constituents were already studied; masks and shared residues
are correlated. This is neither independent-family accuracy qualification nor
selection of pairs after inspecting successful reconstructions.

## Exact execution and independent evidence

The original frozen `research_rosetta_joint.py prepare/score` performs the same
bounded Cartesian search and two-target-view full-pose scoring. `ref2015`, installed
Rosetta revision `1a56185c2592611dec4c9c75ddc9468cd2227c1f`, seed 7, fresh-H imports,
physical/virtual distinction, bounded terminal OXT exceptions, source/topology
checks and complete-feasible-pool selection are unchanged. No packing, relaxation,
reflection, H optimization or post-hoc restoration is introduced.

The original `audit_rosetta_joint.py` replays saved candidates without generating
coordinates. `scripts/audit_rosetta_contacts.py` independently reconstructs the
input-CA ranking and every mask, verifies all unselected records, and supplements
that readback with pool diagnostics, native geometry and cross-generated distances.
[Portable evidence](rosetta-contact-controls-evidence.json) embeds the complete
joint readback and supplemental controls, bound to source/run/code hashes.

| Measure | Count |
|---|---:|
| Enrolled / admitted masked cases | 12 / 12 |
| Source-feasible Cartesian combinations | 26 |
| Raw joint rejections / extra rounded rejections | 3 / 0 |
| Eligible combinations / canonical generated-only `.dat` sidecars | 23 / 23 |
| Native poses, never eligible | 12 |
| Scored full poses / exact target-view imports | 38 / 76 |
| Selected private PDB/`.dat` pairs | 12 |
| Two-varying-axis matrices | 0 |

All 38 scored poses pass both exact physical-heavy/topology imports and stereo
readback. Diagnostic rejected combinations are scored and retained, but cannot
win selection or acquire a materialized `.dat`. Observed records remain unchanged.
All 23 sidecars and all 12 private selected pairs match canonical eligible outputs.
No native/scorer-created H, OXT or virtual atom becomes generated ownership.

## Actual new–new clashes, not a coupled-score success

For THR43/THR51, both constituent proposals pass source-only geometry. Adding the
other target's generated atoms causes these raw rejections:

| Mask | Indices | Conflicting generated atoms | Distance |
|---|---|---|---:|
| Distal/distal | [0, 0] | CG2–CG2 | 2.455 Å |
| Distal/distal | [2, 0] | OG1–CG2 | 2.550 Å |
| Whole/whole | [0, 0] | CG2–CG2 | 2.689 Å |

Both independent contextual audits identify each conflict. The source-feasible
pool remains complete, including these rejected combinations. This supplies real
case evidence for the joint geometry gate beyond synthetic tests. However, the
original **independent per-target clearance winner already avoids all three**;
there are zero independent-winner joint failures. Do not claim a rescue, improved
coverage, necessity of energetic coupling or a general packing advantage.

Every case has an axis with exactly one source-feasible candidate. The enrolled
finite-difference diagnostic consequently labels all 12 matrices **uninformative**,
not proof of absent physical interaction. Across 24 case-target instances there
are 42 source-rejection finding occurrences, all `generated_clash`; correlated
instances/findings are not independent residues or 42 independent failures.

Score changes three raw-joint-clearance choices, without increasing coverage:

| Pair/mask | Previous pooled RMSD | Score-selected pooled RMSD |
|---|---:|---:|
| THR89/VAL92 distal/distal | 1.782 Å | 0.238 Å |
| THR89/VAL92 whole/whole | 1.464 Å | 0.196 Å |
| VAL2/THR40 distal/distal | 1.738 Å | 0.344 Å |

Matched 12-case exact-name medians are **0.233→0.202 Å**, pooled over generated
atoms of both targets. The first two share a pair; the third includes the previously
known VAL2 effect. Native RMSD is reporting only and never selects a candidate.
No template comparator was run here; no superiority claim follows.

## Native geometry versus generator policy

Native coordinates are not staggered-grid proposals. An initial supplemental
probe (`contact-supplement-01.json`) passed them without torsion metadata to the
generator audit and obtained 0/12 accepted materializations. Its retained
`torsion_metadata_or_grid_mismatch` findings are **protocol incompatibility, not
12 native geometry failures**. This did not affect enrollment, candidates, scores
or selection. The corrected reporting-only audit uses the **existing comparator's
measured-χ path**, including both targets' native additions in the context. It
does not alter steric, bond, angle or stereo limits, or permit generator admission.

Under that geometry-only control, **9/12 pass raw and rounded checks**. The three
failures are correlated THR89 masks: deposited THR89 OG1–ASN87 OD1 is 2.399187 Å,
giving clearance −0.640813 Å against the unchanged −0.6 Å cutoff. This is source
geometry, not reconstruction damage; no post-hoc threshold relaxation or label
of a proven chemically impossible contact is justified. All native poses still
pass the separate exact import/stereo checks. Generator materializations remain
0/12 as expected; native controls never become eligible outputs.

This is another reason to keep geometry feasibility, native accuracy, grid policy
and scientific qualification separate. Better RMSD does not calibrate the gate.

## Reproduction, verification and next gate

Retained project-local roots: `tmp/pdbfixer-research/contact-enrollment-01/`,
`contact-plan-01/`, `contact-score-01/`, `contact-audit-01.json` and supplemental
reports. Selected files are private research artifacts, not public bundles or
GUI-registered outputs. Use exclusive new directories for new runs:

```bash
PYTHONPATH=src:. python scripts/research_rosetta_contacts.py --output <enrollment>
PYTHONPATH=src:. python scripts/research_rosetta_joint.py prepare --enrollment <enrollment> --output <plan>
# Existing user-authorized licensed PyRosetta interpreter:
PYTHONPATH=. python scripts/research_rosetta_joint.py score --plan <plan> --output <scores>
PYTHONPATH=src:. python scripts/audit_rosetta_joint.py --plan <plan> --scores <scores> --output <joint-audit.json>
PYTHONPATH=src:. python scripts/audit_rosetta_contacts.py --enrollment <enrollment> --plan <plan> --joint-audit <joint-audit.json> --output <evidence.json>
```

**615 focused tests passed in 26.53 s**; targeted Ruff passed. Tests prohibit
outcome preview during enrollment, prove sidechain-coordinate independence,
reject omitted cases and retain true native clashes while separating grid metadata
from measured geometry. This is not a full repository test or new-family benchmark.

Next gate: separately declare a bounded richer-sampling study on **all frozen
pairs**, not only promising two-axis examples. Establish whether coarse staggered
sampling or the steric policy collapses diversity, without changing observations,
using truth for selection or weakening stereo/clash rules to obtain a positive
coupling result. Calibrated context chemistry remains a separate research need;
an approximate interaction graph and production integration are still unjustified.

Final joint and supplemental readbacks reproduce byte-identically without any
coordinate regeneration. Package mypy, generated CLI/GUI checks, agent maps, old
evidence and all 49 fixture hashes also pass. Final artifact checks are retained
in `tmp/pdbfixer-research/validate-contact-final.py`.
