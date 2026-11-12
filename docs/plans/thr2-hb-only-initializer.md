# Prospective HB-only initializer counterfactual

Continuation of the [frame algebra controls](../research/thr2-initializer-frames.md)
after `1e5d517`, on `research/thr2-hb-only-initializer`. This is a separately
frozen research experiment, not production initialization or a main merge.
Documentation scope: this plan, frozen JSON, complete scalar audit/results and
canonical handoff; validate links/digests and keep coordinates/traces private.

## Question and complete enrollment

Does removing **only THR2 HB's saved perturbation**, before H-only minimization,
preserve a compatible H geometry and alter the ensuing Cβ branch? Global noise
removal is excluded: the preceding mathematical zero-offset control has 139
exact same-parent H coincidences per frame. Do not propose it as a physical arm.

Sources are all three original fresh THR2 captures: unseeded, seed 1, seed 2.
For each, run both **unchanged** and **HB-only-unperturbed** starts once in both
captured-checkpoint OpenCL and Reference XML modes: **12 paths / 36 native
minimizer calls**. Each path executes the original H-only stage and two original
whole-system phases in sequence. These are six deterministic input proposals in
two numerical modes, on one poor-heavy context, not twelve independent samples.

HB-only intervention uses the original heavy-neighbor-derived 0.1 nm direction
from the saved pre-H input, without its random perturbation. Retain all other
314 initialized H and all 327 heavy coordinates exactly. No new atoms, random
draws, donors, orientations, radii, strengths, FFs or platform search. This
changes HB's point, not one isolated force term or torsion alone.

## Source/bridge admission and unchanged physics

- Verify original captured artifacts, runtime provenance and full atom/bond order.
  H/whole residue-label differences may only be canonical parent/variant aliases
  under `ffutils.variants.VARIANT_TO_PARENT`; no hand-rolled variant chemistry.
- H-only: exact original System and 1.0 kJ/mol/nm / 50-iteration call; heavy
  particles pinned and only the original 315 H mobile. Whole phases: original
  Systems and 10 kJ/mol/nm / 1,000 iterations, no solvent. Phase 2 uses the exact
  captured same-anchor divide-ten restraint transition, not reanchored weights.
- Chain stages from their **own returns**, never splice a captured native return
  after a counterfactual. Apply only the exact historical H-repair bridge before
  whole phase 1, not the newly implemented production O–H guard. Saved-only
  admission already reproduced all three original H→whole arrays exactly with
  zero historical repairs; counterfactual repair outcomes remain reported.
- No beta/alpha wall addition, reflection, emergency rescue, optimizer or FF
  retuning. Existing captured native terms remain intact. No candidate PDBs or
  public structural outputs: negative geometry is retained as research evidence.

## Exact masks and numerical transport

Reference requested/applied coordinates must match exactly. OpenCL may quantize
a newly computed point; each coordinate's bound is **half the larger adjacent
float32 spacing plus 1e-15 nm**, calculated before observing its applied value.
All copied initial coordinates and all H-stage heavy coordinates must remain
exactly unchanged. Store requested and actual arrays and measured deltas; no
post-outcome tolerance expansion or platform switch. This does not weaken the
preceding exact-transplant contract: that older study copied already representable
captured coordinates, whereas this intervention computes one new point.

For both numerical modes, verify exact inherited atom/bond identity, source
System bytes and phase continuity. Require unchanged Reference H-only endpoints
to equal their previously retained Reference saved-call replays exactly. Do not
require new chained whole-phase Reference endpoints to equal old replays that
started from native OpenCL H returns. OpenCL endpoint variability is reported.

## Observation, limits and retention

Reporter captures positions/gradients/native statistics without requesting early
exit. Independently re-evaluate physical mobile gradients/objectives at every
saved callback, and report the last-callback/return comparison and native
gradient criterion. Measure heavy Cα/Cβ, RDKit CIP, and the distinct HB
four-ligand volume at each stage; do not conflate volume sign with CIP. Retain all
sign crossings, clashes, bad H geometry, failures and nonexecutions. Neither
rescue nor failure changes enrollment or permits an extra seed.
Maximum O–H and minimum same-parent H–H distances are descriptive node metrics;
complete raw coordinates retain other contacts for inspection. These metrics
are not a new general clash-acceptance gate or a replacement for the old bridge.

Bounds: **4,096 callbacks per phase**, **1,200 s/path**, **7,200 s aggregate**,
at most **36 declared native calls**; no retries or extensions. Preserve partial
traces on failure and kill timed-out process groups. Complete callbacks/native
arrays/XML/checkpoints remain private under
`tmp/pdbfixer-research/thr2-hb-only-run-01/`; tracked outputs contain only hashes,
identities and scalar measurements. Full original sources/failed trajectories
and older negative results remain immutable.

Interpretation is finite-context initializer dependence, not broad qualification,
uniform-noise superiority, an original-lost-run recreation, a general replacement
H generator or production authorization. If this intervention is favorable,
independent contexts and a separate production design/approval gate still remain.
