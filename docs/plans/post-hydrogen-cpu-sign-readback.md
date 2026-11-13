# Prospective source-sign/float32 CPU readback attribution

## Hypothesis and limits

Read-only inspection of exact-revision CommonMinimizeKernel.cpp after the local-
policy native experiment returned identified a separate testable inconsistency.
In `evaluateCpu`, delta is `(positions[y]-positions[x])/r`. For penalty
`0.5*K*(r-r0)^2`, the Cartesian gradient should add `+K*(r-r0)*delta` to y
and subtract it from x. The inspected source subtracts from y and adds to x.
Its `hostGrad` is then uploaded to the single-precision device grad buffer,
which the reporter later downloads. This is a **source hypothesis**, not
independent attestation of installed native code or an accepted physical gradient.
Do not patch OpenMM, retune physics, weaken original tests, or call matching an
objective-inconsistent source formula a scientific pass.

## Fixed enrollment before queries

Use ALL 26 originally failed phase-1 callbacks from the `9a0aeba` secondary-
readback plan, unchanged exact System/Integrator/State/coordinates and integration
force groups. These callbacks all completed both previous secondary readbacks
but passed the complete physical envelopes 0/26 on each backend. Preserve those
failures and every local-policy native outcome.

One default CPU Context, one force/energy query at each of the 26 states, no
native minimizer or H call. Freeze retained exact-revision upstream source bytes,
old plans/evidence and new helper/test hashes. Exactly three predeclared formulas:
(A) physical gradient + analytic positive penalty gradient (original audit);
(B) physical gradient - analytic penalty gradient (source-sign hypothesis);
(C) float32(B) converted back to float64 (source-sign plus upload hypothesis).
Keep unchanged the original per-component absolute/relative envelopes for all
three. No single-precision arithmetic allowance, no fitted offsets or extra
variants. Compare full arrays, not selected convenient components. Quantization
is an explicit prospective hypothesis, NOT a widened tolerance. Energy/penalty
comparisons use the original bounds, and all finite-difference physical tests
must continue to prefer A over B when constraints have nonzero residual.

Record scalar residuals and failed-component counts for ALL states/formulas,
identity-bound witnesses at the largest envelope excess, and whether failures
touch constrained atoms. Retain actual CPU force/penalty/native arrays privately
so the attribution is auditable without re-querying. Mark missing/unavailable
queries and budget nonexecutions explicitly. Maximum 300 seconds / one Context /
26 force queries, no retries, new trajectories, RNG draws or optimizations.

If C matches while A fails, this supports a source-consistent explanation of
reporter gradients; it DOES NOT restore physical-gradient attestation. A wrong
constraint-gradient sign is a substantive native-objective concern, not harmless
rounding. Installed-binary/event causality and upstream reporting require a
separate controlled reproducer/approval scope; no external issue is filed here.
No production capability may be inferred.

Documentation order: this contract + frozen plan, complete scalar readback and
private arrays, then local-policy result/validation and canonical handoff update.
Pending links none; do not publish coordinate/RNG/State arrays. Existing native
contracts, false old gate and failed tolerances remain unchanged.
