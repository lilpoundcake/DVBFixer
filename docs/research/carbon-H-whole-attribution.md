# THR24 snapshot forces: the guard is active, but endpoints are not stationary

Status: **retrospective saved-only attribution**, not an optimizer-trajectory
replay, repair or production qualification. Follows the
[complete whole factorial](carbon-H-whole.md). The original unforced 1CRN THR2
Cβ inversion remains a separate unresolved event.

## Complete matched readback, without new optimization

[Independent evidence](carbon-H-whole-attribution-evidence.json) retains all
**24 THR24 paths / 72 saved states / 216 focal measurements**: both original
seeds, raw/beta-only/full-generated conditioning and all four H arms. Target Cα,
target THR-Cβ and context GLN16 Cα are measured at initial/phase-1/phase-2 states.
The focus is retrospective to already-known near-planar signs; neither negative
endpoints nor favorable controls are selected out of the matched denominator.

Fresh Reference Contexts read the original native System and a measurement-group
clone at identical bound nm coordinates. Captured global parameters must match
native defaults. Restoring group metadata must restore the exact original
serialized System. No minimization, integration, H generation, force ablation,
phase/strength/seed search or coordinate repair occurs.

All original snapshot energies reproduce. Independent all-particle force sums
and signed-volume projections close to the original native System; maximum
absolute component/projection error **1.907e−6 kJ/mol/nm**, maximum energy-sum
error **1.490e−8 kJ/mol**. These are reported absolute errors, not zero. Independent
Cartesian finite-difference gradient error is at most **9.154e−11 nm²**.
All full native/group force arrays remain private and checksummed; no upstream
force-field or empirical rotamer tables are redistributed.

Projection is onto the normalized four-atom Cartesian signed-volume gradient
(CA/N/C/CB, or CB/CA/OG1/CG2), units kJ/mol/nm. Positive means an instantaneous
volume-increasing force component. It is **not** a quasi-Newton step, elapsed-time
velocity, integrated work or proof of a particular unrecorded crossing.

Plan SHA256:
`81a5588e2cc1bac6471841cc31f38b2e9d6b670e57ed3580ea02e0ffbaf78d86`.
Private plan/run/audits: `tmp/pdbfixer-research/carbon-H-whole-attribution-*`.

## The guard exists and agrees with its analytic derivative

All **144 Cα measurements** (target/context over 72 states) have the native guard;
**nine are active**. Every observed guard projection agrees with the independent
analytic derivative. Target Cβ has no Cα guard, as expected. There is no missing
guard, hidden parameter modification or sign/unit mismatch established here.

For raw and beta-only THR24 input, inherited target Cα reference volume is only
**0.025934 Å³**, with `qmin=0.2`; the absolute activation floor is therefore
**0.005186748 Å³**, not the separate diagnostic 0.5 Å³ descriptor. A severely
flattened positive volume just above that floor receives no guard force. In
full-generated conditioning, the target reference is **2.555116 Å³** and the
floor **0.511023 Å³**. Context GLN16 reference is **2.528406 Å³**, floor
**0.505681 Å³**. These are original frozen parameters, not newly proposed settings.

## Negative endpoints already have enormous volume-restoring force

For THR24 / seed 2 / raw, phase-2 force projections are:

| Measurement | Unchanged | Methyl zero | CH₂ | Combined |
|---|---:|---:|---:|---:|
| Cα volume, Å³ | +0.009691 | +0.017269 | **−0.006118** | **−0.006974** |
| Guard active | no | no | **yes** | **yes** |
| Guard projection | 0 | 0 | **+11,554,131.899** | **+12,438,578.149** |
| Native total projection | −270.530 | −351.377 | **+11,553,861.831** | **+12,438,326.631** |

Force units are kJ/mol/nm. The positive guard dominates the negative endpoints,
while the remaining force sum is only approximately −270/−252 along this axis.
All selected whole Systems have **zero constraints and zero frozen particles**.
Thus these large nonzero forces exclude treating the endpoints as converged
stationary minima that inherently prefer D. This does **not** identify whether
the native optimizer stopped at its iteration cap, encountered a line-search
failure or followed some other internal path. No iteration expansion is tried.

The saved phase-1 endpoint in **seed 1 / beta-only / methyl-zero** is also
negative (**−0.002076 Å³**), with restoring guard projection **+7,205,585.536**;
phase 2 ends positive but still flattened (**+0.007662 Å³**). This intermediate
finding was already present in complete whole evidence and must not be erased by
its positive final sign. Across all 408 whole states there are exactly those
three negative-Cα snapshots and no negative THR-Cβ snapshots.

## Other components and limits

At seed-2/raw phase-2 endpoints, native angle-force projections are approximately
**−1109/−1163/−1130/−1135**, torsion **+817/+840/+829/+850**, and nonbonded
**−170/−172/−171/−173 kJ/mol/nm**, respectively. Positional restraints contribute
approximately +159–163; bonds and CMAP supply the remaining smaller terms.
These are whole force-class projections, not isolated pair causes or evidence of
which line-search direction was actually taken.

Initial THR24 raw Cα projections are positive in all four arms, while phase-1
endpoints vary around small positive/negative totals. The active phase-2 guard
force therefore cannot be used to invent an unobserved monotonic force history.
The GLN16 active guards and all full-generated matched contexts remain reported.
Inherited poor-heavy geometry, finite relative guards and fixed-budget optimizer
behavior must remain distinct from original THR2 causality or production safety.

## Validation and next boundary

Nine new focused attribution/audit tests and seven existing signed-volume/guard
tests pass, including finite differences at positive/zero/negative near-planar
volumes, exact original native readback with minimization forbidden, source
energy/global/identity refusals and immutable native System serialization.
Targeted Ruff passes; independent evidence publications agree byte-for-byte.
The full repository fast suite through this attribution passes **3,055 tests**,
three optional Torch/Biotite skips and 33 slow deselected. Retained log:
`tmp/pdbfixer-research/carbon-H-whole-attribution-fast-suite-01.log`.

Next: freeze a separate **observational replay** of these same 24 paths with
native whole-minimizer callbacks, requiring exact original starts and both
endpoints. Preserve the original two Reference/10/1000 phases and phase-2 weight
reduction, no new seed, force, H regeneration, reflection, retry or budget tuning.
Record callback coordinates, gradients and native objective metadata without
requesting an early exit. These would be computational replays, not new samples.
Until exact trajectory transport is established, do not claim a termination or
crossing mechanism from these sparse snapshots.

Code: `scripts/attribute_carbon_H_whole.py`,
`scripts/audit_carbon_H_whole_attribution.py`.
