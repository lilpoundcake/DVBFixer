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

## Native outcome: local projection is not sufficient

All three paths executed once and returned both phases: **six minimizer calls,
12,128 callbacks**. The [complete all-node native audit](post-hydrogen-local-noop-native-audit.json)
has no audit errors, no missing phases and no replacement paths. Phase Systems
are byte-identical to the captured source; requested/applied transport and own-
return continuity pass. Both source states remain enrolled. All phases retain
DAN/EPE and 213 constraints; no new solvent, H placement, repair, reflection,
post-hoc projection or installed/native force modification occurred.

| Path | Callbacks phase 1 / 2 | Final incompatible moved HA | Phase-1 negative-Cα nodes | Physical-envelope failures phase 1 / 2 |
|---|---:|---|---:|---:|
| baseline projection | 1960 / 1571 | **A:286/HA** | 7 | 26 / 0 |
| guarded unmodified control | 2312 / 2329 | A:0/HA, A:283/HA, A:286/HA | 8 | 26 / 0 |
| guarded projection | 2291 / 1665 | **A:286/HA** | 9 | 26 / 0 |

The two projection paths improve final SER0/HA and HA283 relative to the
unmodified guarded control (the reused original baseline also retained a final
SER0 conflict). Nevertheless **0/2 proposal paths pass the prospective geometric
criterion**. HA286 crosses back to an incompatible sign at phase-1 callback31
in both, and is incompatible at both phase returns. Its final heavy normalized
volume is still positive but small (~0.0360/~0.0369), while its HA four-ligand
volume is strongly negative (~-0.9773/~-0.9861). Normal final C-H lengths of
~1.0852/~1.0857 Å do not rescue this wrong-side geometry.

SER0 has four incompatible phase-1 nodes in each projection path, first at
callback26 when its heavy sign is temporarily negative. All three paths contain
transient D-Cα; every final returned state has **zero** D-Cα. No phase-2 negative
Cα node appears. These are retained private downstream States, not published
production structures, and no unplanned emergency repair was applied to hide
the prospective all-node failure.

PRO379/HA stays compatible at every recorded node and relaxes from its warned
1.650683 Å start to ~1.0909 Å in both guarded paths. This shows that the compatible
no-op policy did avoid the old admission blocker locally. It does not establish
comprehensive geometry correctness or justify accepting every elongated bond.

All starts remain astronomically high-energy (~2.3454e13 kJ/mol). Under unchanged
original physics envelopes, **each phase 1 fails at 26 callbacks; all phase-2
callbacks pass**. The new six-phase physical qualification therefore fails
separately from the HA286/all-node geometry failure. The source-sign study above
was confined to the ORIGINAL 26 callbacks; do not claim it directly attested the
78 new failed callbacks or proved that CPU-sign behavior caused HA286 inversion.
Constraint projection/return motion remains explicit and constrained stationarity
uncertified. No force/threshold/acceptance endpoint was retuned after these data.

## Validation and decision boundary

Final full fast suite: **3,292 passed**, three optional Torch/Biotite skips,
33 slow deselections. All **49 focused controls** (23 new in this branch),
Ruff, strict-target mypy, generated CLI/GUI specs and agent maps pass. The
[validation/quality record](post-hydrogen-local-noop-validation.json) separates
these engineering passes from **failed native geometry and physical qualification**.
Variant 1/2, core strict proposer, old API adapter, production and fixtures are
unchanged. Production/default registration and main merge remain disabled.

The local compatible-noop policy is a valid bounded atomic proposal on these
two retained states, but **not a sufficient post-H safeguard through whole
relaxation**. Do not enable it publicly, retry this cohort, silently add a
post-whole projection or change native constraints/forces. The next work needs
a separately chosen research scope: an isolated native constraint-sign
reproducer/dependency investigation, revisiting the earlier initializer route,
or a separately specified later-stage geometry policy. No such new native
experiment, dependency patch or public issue is executed in this record.
