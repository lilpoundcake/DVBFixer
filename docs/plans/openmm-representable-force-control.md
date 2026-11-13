# Prospective representable-force fallback threshold control

## Why a new series, and what stays negative

First synthetic freeze `287eacb` completed all12 calls /162 callbacks but its
primary signal and all-controls gate are FALSE. Excessive k=1e14 stress gives
incorrect OpenCL stress gradients even without constraints, rather than the
intended CPU-sign witness. One low-stress callback also fails the old physical
energy bound. Preserve the complete old result and source, not a tuned replacement.

The pinned `minimize.cc::convertForces` defines `limit=0x4000000000000000`
and force scale2^-32, so it requests CPU fallback above **2^30** force units,
BEFORE the signed fixed-point capacity2^31. The original series used the latter
as an observation criterion and an excessively out-of-range stress. Inspecting
source and those negative observations motivates a DIFFERENT, prospectively
frozen trigger within representable range, not retrying unchanged cases until
green. Source event is still not directly observed. No installed dependency/
production change or broader structural readiness is authorized.

## Fixed additional enrollment

Reuse the same three-particle bond/constraint mathematical model and derivative
audit, but particle2 starts at (10,10,0) nm; new stress levels1e7 (below trigger)
and5e7 (representable trigger). Initial maximal physical force is3e8 or1.5e9,
respectively, both below signed capacity2^31; the latter is above detector2^30.
With a unit-scale first optimizer step the trigger should remain above threshold
at the first callback. Verify this from the data, not as an assumption of success.
Original pair separation0.100005 nm, bond k10000/r_eq0.12, optional constraint0.1,
positive12-Da masses, same Verlet/no-integration/tolerances/platform properties.

Twelve NEW cases: OpenCL/CPU/Reference × constraint on/off × below/trigger.
Exactly one minimizer call/case with **maxIterations1**, tolerance10, native
constraint tolerance1e-5. The one-iteration setting is explicitly a diagnostic,
not convergence/physical qualification; preserve any native outer constraint
iterations. At most256 callbacks/case,30 s/case,360 s total. This adds12 synthetic
calls to the prior12, with no retries, additional constants/initial points,
platform substitutions, CPU-sign patches or tolerance widening. Stop after
accounting these12 even if the trigger is inconclusive. All failures/nonexecutions
must remain represented. Separate process/cache/temp ownership is unchanged.

## Predeclared criteria and comparisons

Use the unchanged original analytic/termwise-FD and energy/gradient bounds, with
both strict and predeclared OpenCL constraint-arithmetic allowance fields.
Three hypotheses remain physical+penalty, physical-penalty and float32 of the
source sign, each evaluated for all components. No fitted offsets or precision
allowances beyond the original contract.

Primary witness must be a constrained OpenCL-trigger callback whose reported
largest gradient lies strictly between2^30 and2^31, nonzero penalty gradient,
and physically correct finite differences/physical and penalty energy envelopes.
Physical+penalty must FAIL even with the original arithmetic allowance; at least
one source-sign hypothesis must PASS the original strict component envelope
(not merely a new allowance). All other11 controls must have callbacks, completed
calls and physically correct original-envelope gradients/energies/finite
differences. Preserve every callback rather than cherry-pick a successful node.

Matching source sign is an installed PUBLIC-REPORTER objective inconsistency,
not proof of binary/source identity or direct internal CPU-fallback event.
Describe any new result as local synthetic behavior only. This does not
rehabilitate structural HA286 geometry, original physics attestation, the first
negative toy series, or establish a corrected implementation. A dependency patch
or external issue remains a separate user decision.

Documentation: this separate amendment/plan/test freeze, complete second-series
scalar evidence, updated first-series result note, final validation/handoff.
Pending links none. Bound exact kernel/source bytes and preserve first hashes.
Raw toy coordinates, Systems/States/checkpoints remain in ignored project tmp;
no protected external native-library discovery is retried or bypassed.
