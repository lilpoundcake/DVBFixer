# Downstream methyl control: limited heavy-gate gains, fixed-phase regressions retained

Status: **one-protein native-phase mechanism control**, not complete complex
acceptance or a production repair. Follows the
[tetrahedral H-stage experiment](tetra-methyl-control.md). No production code,
force/threshold/default, empirical table or public completion capability changed.

## Faithful post-H bridge and all 72 native continuations

[Complete evidence](methyl-whole-control-evidence.json) retains **72/72 completed
cells**, 24 original contexts × unchanged / fixed-phase tetrahedral /
initial-energy-selected tetrahedral H start. All 216 initial/phase-1/phase-2 states
remain reported. No retries, execution blockers, timeouts or budget nonexecutions.
These remain four masks, two H seeds and three source-geometry arms on **1SHG**,
not a broad or independent-protein qualification.

Do not copy the H-only endpoint straight into a full-system interpretation:
the original production path applies `repair_misplaced_hydrogens` between them.
`research_h_stage_transport.py` reconstructs the original atom/bond topology
without renaming IDs and applies that **existing helper exactly once**. A separate
preflight reproduces all **24 original whole-system starts exactly in native nm
coordinates**, including H, not merely PDB rounding or Å-level tolerance. Heavy
coordinates remain exact through this bridge. All grid/intermediate artifacts
from the prior H-stage experiment stay bound and immutable.

The continuation uses the original saved full Systems, State and integrator,
**Reference / no solvent / tolerance 10 kJ/mol/nm / 1,000 iterations per phase**.
Only the original positional weights are divided by ten for phase 2; all original
anchors/reference volumes and other force parameters remain unchanged. Both
saved endpoints of all **24 unchanged controls reproduce exactly**, all heavy
atoms and H. Each published state independently checks its serialized System
against the original corresponding phase.

Plan: `tmp/pdbfixer-research/methyl-whole-plan-01.json`, SHA256
`1ff54b50dd4f879dedccb94620a626a4c53645ecf7cfcd4e582af3d2dd0b24a9`.
Private coordinates/Systems/logs are under `methyl-whole-run-01/` beside it.
Budgets were 1,800 seconds per cell and 7,200 seconds total. No PDB, accepted
candidate bundle or production result is published; private paths do not transfer
through Git.

## What this gate does and does not measure

Reuse the unchanged actual-frame raw and PDB-rounded **target heavy geometry +
global Cα/THR-Cβ stereo/shape** research gate, with no outcome-derived thresholds.
It is not an exhaustive all-heavy/all-H clash, chemistry, energy-convergence or
biological-assembly certificate. All H-bond and H-containing harmonic-angle
deviation summaries are separately retained, not silently folded into a new gate.
There is no held-out native-coordinate/RMSD comparison or observation restoration.

## Limited overall gains

Raw and rounded pass counts coincide:

| Original source-geometry arm | Unchanged | Tetrahedral zero | Tetrahedral selected |
|---|---:|---:|---:|
| Raw | 2/8 | 2/8 | 2/8 |
| Distal-only (`beta_only`) | 3/8 | 3/8 | 4/8 |
| Full generated-heavy rebuild | 5/8 | 6/8 | 6/8 |
| **All** | **10/24** | **11/24** | **12/24** |

Both canonical arms prevent the retained **observed GLN16** failure for
full-generated THR24 seed 1. The selected arm additionally changes distal-only
THR24 seed 2 from target-Cα failure to a pass (final target volume
0.009561→2.480606 Å³). This does not cure every target: the selected arm still has
**six THR Cα/angle failures and six VAL23 branch-name orientation failures**.
VAL's equivalent methyl-name convention failure is not a chemical Cβ inversion.

For the focal full-generated THR24 seed-1 case:

| Final GLN16 measurement | Unchanged | Tetrahedral zero | Tetrahedral selected |
|---|---:|---:|---:|
| Heavy Cα volume | 0.482378 Å³ | 2.544202 Å³ | 2.545041 Å³ |
| N–CA–HA | 85.944° | 109.453° | 109.432° |
| HA–CA–CB | 85.697° | 107.914° | 108.001° |
| HA–CA–C | 89.067° | 108.270° | 108.224° |
| CA–HA radius | 1.087843 Å | 1.086718 Å | 1.086862 Å |
| Global maximum H-angle deviation | 23.803° | 6.202° | 6.194° |
| Original 469-heavy-atom RMS movement | 0.219522 Å | 0.219133 Å | 0.218328 Å |
| Original heavy maximum movement | 1.077288 Å | 1.098176 Å | 1.084132 Å |

The later native relaxation resolves the focal H-stage bond compression/angle
strain as well as its heavy-Cα flattening. But original observations **do move**;
this is not an observation-immutable coordinate repair. Only the preceding H
initialization/transport preserves heavy coordinates exactly.

## Fixed-phase regressions: aggregate pass counts are insufficient

The fixed-phase arm has two important adverse findings in already rejected cases:

1. **THR24 / distal-only / seed 2:** new **observed TYR15 Cα near-planarity**,
   signed volume **+0.425081 Å³**. The target THR24 still fails, so this new context
   failure is invisible in the overall 11/24 pass count. It is positive-volume
   shape damage, not D-Cα.
2. **THR24 / distal-only / seed 1:** target Cα crosses to **negative volume
   −0.014602 Å³**, whereas its unchanged counterpart remains near-planar positive.
   This is a retained **D-Cα native-phase finding**, not an accepted repaired
   candidate. Zero such final negative-Cα cases occur in unchanged/selected arms;
   fixed-phase has one. No final negative THR-Cβ was found in any arm.

These are deliberately the original **two native phases before production's
post-phase-2 force-reflect safety net**. No emergency reflection was added to hide
a research failure. Production still has its unchanged zero-D-Cα output safeguard;
this study neither removes nor disproves that final safeguard. No invalid PDB is
published. The fixed-phase construction must not be called a safe repair merely
because its aggregate gate count improves or its H-only endpoint avoids inward HA.

The selected arm has no newly appearing global stereo/shape findings in this
cohort, but 24 correlated starts on one protein cannot establish general safety.
It must not inherit a production capability, default or broad qualification.

## Validation and continuation boundary

Ten focused transport/native-phase/H-angle tests pass; targeted Ruff and agent-map
and generated CLI/GUI checks pass. The full fast suite through this continuation
passes **2,906 tests**, three optional Torch/Biotite skips and 33 slow deselected.
Independent publications regenerate byte-identically. Original source/failed
outputs and earlier numerical evidence remain immutable.

Next: separately enroll a source-diverse, graph-defined methyl initialization
control with explicit generated-H ownership and fixed input/energy policies before
seeing its outcomes. Retain failed admissions and adverse contexts; do not tune
phases/seeds to the known positive case or enable production. Original unforced
1CRN THR2 Cβ causality and the empirical-table production/license gate stay open.

Code: `scripts/research_h_stage_transport.py`,
`scripts/research_methyl_whole_control.py`, `scripts/audit_methyl_whole_control.py`.
