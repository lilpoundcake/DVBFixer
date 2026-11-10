# A one-HA start change prevents the retained GLN16 flattening endpoint

Status: **case-specific controlled initial-condition evidence**, not a production
repair or an explanation of the original 1CRN THR2 Cβ inversion. Follows the
[saved-state angle attribution](stereo-angle-terms.md). No production code,
backend, threshold, force parameter or public completion capability changed.

## Frozen paired experiment

All 24 retained starts were enrolled, each with `unchanged` and `outward_HA`
replays: **48 declared cells, all 48 executed successfully**, no retries,
timeouts or budget nonexecutions. This is four masks on **one protein, 1SHG**,
two retained H-placement seeds and three source-geometry arms, not 24 independent
proteins. Each matched pair uses the same original saved System, integrator,
State, topology and two-phase nonperiodic **Reference / no-solvent** protocol:
tolerance 10 kJ/mol/nm, at most 1,000 iterations per phase. Only the original
positional restraint weights are divided by ten for phase 2; anchors and
chirality reference volumes remain unchanged. Budgets were 1,800 seconds per
cell and 7,200 seconds for the complete cohort.

The intervention places only **generated GLN16 HA** opposite the sum of the
three unit CA-to-N/C/CB heavy-bond vectors, keeping its initial CA–HA distance.
It is applied to every start, not selected using initial defects or native error.
There is no H regeneration, new force, strength/seed grid, reflection or source
coordinate restoration. Source records are checked to ensure HA was not an
original observation; the captured initial State must match the starting atom
coordinates exactly. Every other starting coordinate remains exact.

The final plan is retained privately at
`tmp/pdbfixer-research/ha-start-plan-01.json`, SHA256
`ec5500f041906f721a5ad10aecebe6c230304f75905687060cfb66ad8346d1ae`.
The first enrollment attempt failed before writing a plan or executing a cell:
the old trace has no `settings` field. The driver was corrected **before** the
plan was frozen to inspect its actual phase/platform records. The explicit
settings match the hash-bound historical CLI default and OpenMM API default;
no scientific outcome or frozen plan was rewritten.

## Independent audit and unchanged replay

[Complete evidence](ha-start-control-evidence.json) retains every cell, all three
saved frames, target/global raw and PDB-rounded geometry, GLN16 HA angles,
original-observation motion and all 24 pair comparisons. The saved-only auditor
was written while sampling ran, before replay results were inspected, using the
predeclared unchanged geometry gates. It independently checks complete atom
identities, the one-H-only intervention, negative-unit-sum direction and preserved
CA–HA radius. Each saved System must serialize identically to the original
phase's System, checking all forces, parameters, anchors and reference volumes.

**All 24 unchanged controls reproduce both original saved phase endpoints
exactly**, for all heavy atoms and hydrogens: maximum difference **0.0 Å**.
Thus the intervention can be compared against a reproduced baseline rather than
an approximate restart. Every paired start preserves all non-HA coordinates.
Raw states, forces and execution logs remain private under
`tmp/pdbfixer-research/ha-start-run-01/`; these paths do not transfer through Git.
The independent audit/publication regenerates byte-identically without replaying.

## GLN16: the retained failing trajectory changes

THR24 / full-generated / source seed 1:

| Measurement | Unchanged | Outward HA |
|---|---:|---:|
| Initial GLN16 Cα volume | 2.528406 Å³ | 2.528406 Å³ |
| Initial HA movement | 0 Å | 2.149694 Å |
| Final GLN16 Cα volume | **0.482378 Å³** | **2.546770 Å³** |
| Final N–CA–HA | 85.944° | 109.488° |
| Final HA–CA–CB | 85.697° | 107.942° |
| Final HA–CA–C | 89.067° | 108.227° |
| Final raw / rounded heavy geometry | fail / fail | pass / pass |
| Original 469-heavy-atom RMS motion | 0.219522 Å | 0.218545 Å |
| Original heavy-atom maximum motion | 1.077288 Å | 1.097194 Å |

For this retained start, changing only HA under the same reproduced physics
**prevents the GLN16 near-planar endpoint**. Together with the independent
H-containing-angle force attribution, this supports an H-start-sensitive basin
mechanism for this particular GLN16 failure. It does not establish a unique
microscopic optimizer path, show that every inward H causes flattening or explain
why the upstream H-placement procedure produced that initial HA.

Minimization does move original observed heavy atoms; the result is **not** an
observation-immutable coordinate repair. The initial comparison to deposited PDB
decimal coordinates also reports floating conversion noise (up to 5.3×10⁻¹⁵ Å
in this case), rather than silently calling that exact source-byte preservation.
The exactness claim concerns the matched replay starts and unchanged endpoints.

## Whole cohort: limited improvement, failures retained

| Source geometry arm | Unchanged raw / rounded passes | Outward-HA raw / rounded passes |
|---|---:|---:|
| Raw | 2/8 / 2/8 | 2/8 / 2/8 |
| Distal-only (`beta_only`) | 3/8 / 3/8 | 3/8 / 3/8 |
| Full generated-heavy rebuild | 5/8 / 5/8 | 6/8 / 6/8 |
| **All** | **10/24 / 10/24** | **11/24 / 11/24** |

Only the GLN16 pair above changes final pass/fail status. All outward-HA final
GLN16 volumes are above the unchanged near-planarity bound (range 2.305469–
2.573437 Å³); the seven other THR Cα geometry failures and six VAL23 branch-name
orientation failures remain, yielding **13/24 outward failures**. This is not a
general THR/VAL reconstruction solution. No native-reference RMSD was computed,
and earlier accuracy and source-refusal evidence is unchanged.

## Validation and next boundary

The driver/auditor have **18 focused tests**, including independent synthetic
one-H geometry, source-observation ownership, exact State checks, undeclared
motion and missing/duplicate/nonfinite identity refusals. The combined new
angle/HA set passes 40 tests; targeted Ruff passes. The full fast suite passes
**2,830 tests**, with three optional Torch/Biotite skips and 33 slow tests
deselected. Agent-map and generated CLI/GUI checks pass. This remains Reference-only research,
not a slow-suite, platform, protein-diversity or production qualification.

The subsequent [original H-placement trace](ha-placement-trace.md) now reproduces
all 24 original H endpoints exactly: the known HA starts all-obtuse and becomes
all-acute during native H-only minimization, with existing heavy atoms immobile.
This localizes the upstream stage; it does not justify installing this one-H
construction as a fix. Keep that provenance question distinct from the
still-unreproduced original THR2 event and the remaining THR Cα failures. The
licensed rotamer-table production gate and disabled public repair capability
remain unchanged.

Code: `scripts/research_ha_start_control.py`, `scripts/audit_ha_start_control.py`.
