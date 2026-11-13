# Private compatible-noop policy and retained-state downstream controls

## Scope and prospective contracts

After the explicit local-policy proposal at `5284081`, the user instructed
continued private work. Work remains on `research/post-hydrogen-local-noop`;
production, installed OpenMM, main, original assertions and all older strict/
variant-1/variant-2 evidence remain unchanged. This is a separate policy, not a
relaxed pass for the [failed public-API experiment](post-hydrogen-public-api.md).

The [offline contract](../plans/post-hydrogen-local-noop.md) and
[two-state frozen plan](post-hydrogen-local-noop-plan.json) were committed at
**`666c745`** before evaluating either scientific state. Documentation order is
contract/plan, outcomes/result, validation and canonical handoff; no unresolved
future links. The native study has a separate contract/freeze below.

## What changed, and what did not

All original new H retain their original ownership. After global identity/graph/
coordinate validation, supported canonical carbon-H centers must have
nondegenerate heavy/four-ligand signs under the unchanged normalized 1e-4 bound.
Compatible centers are exact no-ops. If their distances are outside the old
admission intervals, they receive a structured identity-bound warning, NOT a
claim that their bond lengths or overall geometry are acceptable.

Only sign-incompatible supported new H enter an explicitly separate intervention
mask. The unchanged strict proposer receives that mask over the **complete**
original geometry. Every moved H keeps the original radius, heavy-bond bounds,
formula/sign checks and final full-system nonparent-contact check; untouched
atoms remain present as possible contact partners. No atom truncation, radial
fix, H-only relaxation, new force or tolerance increase is involved. Proposal
and warnings are returned only after all interventions admit atomically.

Original global strict admission is still a different, more conservative policy:
it refuses the compatible but elongated PRO379/HA even though no projection of
that H is needed. The new policy explicitly does not fix that elongated bond.
Production/default registration and the live public-API adapter remain disabled.

## Offline outcomes

The [complete two-state results](post-hydrogen-local-noop-evidence.json) used the
already retained ordinary-repair coordinates, not a repeated repair or new H
placement. Both states proposed and passed a separate NumPy audit. Each contains
2,901 original new H, of which 448 are supported and 2,453 outside this selective
carbon-H rule. These are inventory counts, NOT independent samples.

| Retained H-return state | Projected H | Compatible distance warnings | Audit |
|---|---|---|---|
| baseline | HA at A:0, A:283, A:286 | none | pass |
| guarded | HA at A:0, A:283, A:286 | A/PRO379/HA, **1.650683 Å** | pass |

All observations/no-op coordinates remain exact, Python/NumPy RNG unchanged.
Maximum C-H radius differences are ~2.22e-16/~1.94e-16 nm; minimum final
nonparent distances for moved H are ~0.08448/~0.08617 nm (old minimum 0.05 nm).
SER0/HA becomes sign-compatible with unchanged length 1.13303/1.12309 Å. The
warning at PRO379 is retained; the original strict after-repair gate is STILL
false. Neither source state nor any old result was rewritten.

Nineteen synthetic controls cover exact old/new equality on admitted inputs,
compatible elongated no-op warnings, strict incompatible-distance refusals,
degenerate signs, full-inventory contact rejection, immutable observed H,
canonical aliases, deterministic/RNG-neutral and frame-equivariant behavior.
The independent audit rejects corrupted observation, projection, warning and
intervention-mask outputs. The combined new/existing focused set passed 43 tests.

## Separately frozen native downstream control

The [native contract](../plans/post-hydrogen-local-noop-native.md) and
[three-path plan](post-hydrogen-local-noop-native-plan.json) were frozen at
**`8850683`** before execution. All paths use the exact captured baseline whole-
System physics, including DAN/EPE, water constraints, native restraint anchors
and the exact same-anchor divide-ten phase transition. Admission verifies atom
identity/element/canonical residue chemistry, graph and bond metadata, raw API
selected variants and heavy-coordinate/anchor transport for both source states.

This is explicitly **shared-physics counterfactual research**, not proof that an
independently built guarded ForceField would match or that a live API integration
works. The two source H placements differ; they are never presumed identical.

Three new paths are enrolled, without retries: baseline projection, previously
unexecuted guarded-unmodified control, guarded projection. Each has two ordinary
whole phases; the original baseline-unmodified native trajectory is reused,
not repeated until green. No new H samples, ordinary repairs, extra H relaxation,
reflection, platform substitution, or added solvent. Source physics derives from
the original `minimize --no-solvent --keep-water` run. Native Systems, trajectories
and RNG/State/checkpoint payloads stay private under the project tmp tree.

The original constraint-aware physics audit and fixed component-wise envelopes
remain mandatory. Geometric preservation and physical numerical qualification
are separate outcomes; the old 26 failed gradient callbacks stay failed even
if the new geometric endpoints improve. Constrained stationarity and public
readiness are not inferred from a returned minimizer or a zero-D-Cα structure.

## Separate source-sign attribution of old CPU readback failures

Exact-revision source inspection identified a checkable sign discrepancy in
`CommonMinimizeKernel::evaluateCpu`. With delta `(positions[y]-positions[x])/r`,
the positive penalty derivative adds `+K*(r-r0)*delta` to y and subtracts it
from x. The inspected CPU fallback code does the reverse. The device upload can
also quantize the resulting gradient to float32. This inspection does not attest
the installed binary and is NOT a reason to redefine the physical audit.

The [separate prospective contract](../plans/post-hydrogen-cpu-sign-readback.md)
and [fixed 26-state plan](post-hydrogen-cpu-sign-plan.json) were frozen at
**`5ef747c`** before force queries. The
[complete scalar results](post-hydrogen-cpu-sign-evidence.json) used one default
CPU Context and one energy/force query for each of the original 26 callbacks,
without minimization, H placement, repairs, fitted offsets or tolerance changes.
Both physical and source formulas were declared in advance and evaluated over
the complete gradient arrays with the ORIGINAL component envelopes.

| Formula | Full gradient-envelope passes |
|---|---:|
| physical gradient + positive penalty derivative | **0/26** |
| physical gradient - positive penalty derivative | **26/26** |
| float32 transport of the source-sign formula | **26/26** |

All **16,120** failed physical-gradient components belong to constrained atoms;
no unconstrained components fail that comparison. Source-sign agreement already
passes without requiring float32 transport, so do not claim quantization alone
explains or repairs these failures. Maximum absolute residuals at extreme-force
components (~2.68e5, ~2.62e5 after float32) remain explicitly recorded; passing
here means the unchanged absolute/relative envelope, not bit equality. For
example, the first state's A/HOH439/H2 x component is ~218.096685 for the physical
formula versus reported ~217.583649, outside its ~0.005362 original bound.

Two synthetic finite-difference/comparison tests preserve the physical meaning:
the positive derivative matches energy finite differences; the reversed source
sign does NOT. A native-source explanation is therefore **not physical-gradient
qualification**. This supports a source-consistent CPU-fallback constraint-sign
concern rather than harmless GPU arithmetic alone. Installed binary/event
causality, a controlled upstream reproducer and any external issue or dependency
patch remain separate gates. No external issue was filed and no OpenMM code
was changed. The old 0/26 physical outcomes and strict admission refusal stand.
