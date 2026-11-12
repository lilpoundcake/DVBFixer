# THR2: local H mechanism and bounded matched-start transplants

## Scope and saved-only localization

Continues the [fresh native capture](thr2-native-capture.md) after the user's
instruction to continue research autonomously. Production code, force fields,
minimization settings and safeguards remain unchanged. Documentation scope is
this protocol/results note, a link from the capture note and the canonical
[handoff](../plans/pdbfixer-reconstruction-handoff.md); validate evidence and links
before publication. Original capture enrollment is finished, not extended.

`audit_thr2_angle_terms.py` analytically evaluates every native harmonic angle at
**all three retained whole-phase-1 starts**, with independent Reference force
readback. It does not optimize, place H or mutate forces (measurement groups only).
Finite-difference tests validate the analytic Cartesian derivative. Multi-H
methyl angles remain joint terms, never attributed twice.

| Start | Heavy-only angle Cβ projection | Target HB angle contribution | All angle terms |
|---|---:|---:|---:|
| unseeded | +909.351 | **−2,728.961** | −1,414.324 |
| seed 1 | +909.351 | **+2,978.857** | +4,469.238 |
| seed 2 | +909.351 | **+2,980.424** | +4,456.165 |

Units: kJ/mol/nm projected along the increasing-Cβ-volume Cartesian unit gradient.
Maximum analytic/native Cartesian-force discrepancy is **4.1×10⁻¹² kJ/mol/nm**.
The [complete angle audit](thr2-angle-terms-evidence.json) includes every focal
angle and H grouping. HB is a strong initial-force candidate, not yet proof of
whole-trajectory causality. Other H terms also differ.

## Prospective H-swap protocol

The [frozen plan](thr2-hswap-plan.json) is a byte-identical public metadata copy
of `tmp/pdbfixer-research/thr2-hswap-plan-01.json`; native payloads remain private.

The HB mask is chosen **after saved-only localization, before any transplant
minimization**. This is mechanistic follow-up on one poor-heavy THR2 context, not
an independently selected validation cohort. No new seed or native H placement.

Sources: the three retained original whole starts (unseeded, seed 1, seed 2),
with complete matching atom/bond identity, byte-identical Systems and identical
327 heavy coordinates. Transplant coordinates by exact atom identity only:

- **HB:** just `A/2/(empty insertion code)/THR/HB` (one H).
- **other-H:** every H except that HB (314 H); complementary partition.
- **unchanged:** each of the three starts without modification.

Enroll **all six ordered donor→recipient pairs × both masks**, plus three
unchanged baselines: **15 paths**, each once on captured-checkpoint OpenCL and
Reference XML routing, **30 paths / 60 native full-phase calls** total. No
outcome filtering, additional masks, donors, seeds, repeats or platform search.
Full-H substitution closes exactly to the donor start and needs no extra run.

Both native phases use **10 kJ/mol/nm / 1,000 iterations**, no solvent. Phase 2
chains from that path's own phase-1 return, not the original native phase-2
coordinates. Apply only the exact captured native restraint divide-ten change,
keeping anchors fixed; require serialized phase-2 System equality. No force
addition, H repair, reflection, coordinate restoration or production output.

Every actual starting coordinate must equal its frozen transplant exactly; do
not widen a transport tolerance. Retain starts/callbacks/returns privately.
Reporter copies native coordinates/gradients/stats and never requests early exit.
All outcomes, including D geometry, clashes, failures and nonexecutions, remain
reported rather than repaired or silently discarded.

Bounds: **4,096 callbacks per phase**, **1,800 s per path**, **14,400 s total**.
No retries or extensions. Independently audit exact masks/phase continuity,
System equality, stereochemistry and physical gradient/objective readback.
Unchanged Reference phase-1 controls must match the earlier saved-start replay
exactly; OpenCL endpoint variability is measured, not hidden.

Interpretation: successful bidirectional HB flips would support HB dependence in
these retained matched starts, not a universal mechanism, seed success rate,
production initializer or recovered November 6 trajectory. Complementary masks
test whether that interpretation survives the rest of H geometry. No outcome is
a production-qualification gate.
