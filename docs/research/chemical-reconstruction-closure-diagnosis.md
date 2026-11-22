# CD-only ring-closure diagnosis: WP2Q

Date: 2026-11-22. Baseline `6810f76`, persistent
`feature/chemical-reconstruction`. This is constructor-free retrospective analysis
of the frozen [WP2P comparison](../plans/chemical-reconstruction-comparison-validation.json),
not a repair, new reconstruction trial or a change to its results.

The [protocol](../plans/chemical-reconstruction-closure-diagnosis-v1.md) and
[nine-artifact ledger](../plans/chemical-reconstruction-closure-diagnosis-inputs-v1.json)
were committed as `ab3fa5d` before diagnostic measurements. The user approved
diagnosing the two retained original-proposal rejections. Version stays 0.10.2;
no production/default/public/dependency or data-license gate changes.

## Result: a local chain frame does not solve the other ring edge

Both rejected CD-only proposals use exactly the saved **CA–CB–CG -> CD** chain
placement. Its length CG–CD, angle CB–CG–CD and signed torsion CA–CB–CG–CD are
reproduced in raw coordinates to floating-point precision. Maximum raw absolute
residuals across the three CD-only cells are below 2e-15 A and 6e-14 degrees.
There is no evidence here of an `extend` sign error, atom-ID mismatch or rounding
being the cause.

The other bond to the missing atom is **N–CD**. It is an authoritative graph edge
on a cycle (detected by connectivity after removing that edge), but is NOT the
length edge of the chosen chain placement. N is not one of that placement's three
anchors. The constructor does not jointly solve both incident bond constraints;
the independent validator correctly continues to enforce both.

Pinned authority N–CD length: **1.4866485798600826 A**, tolerance **0.10 A**;
allowed interval **[1.3866485798600825, 1.5866485798600827] A**.

| CD-only source | Raw N–CD (A) | Rounded N–CD (A) | Raw undershoot below lower bound (A) | Original own audit |
|---|---:|---:|---:|---|
| 1TIM A/PRO2 | 1.377740014752 | 1.378295686709 | 0.008908565108 | rejected |
| 2CI2 I/PRO25 | 1.354268366481 | 1.354521686796 | 0.032380213379 | rejected |
| 1CRN A/PRO5 | 1.536235114242 | 1.536006184883 | none | passed |

For BOTH rejected sources, the **only original-own failure is
`authority-bond-geometry` on N–CD**, in both raw and rounded phases. All measured
fixed-only original-own terms pass. The violation predates PDB rounding; rounding
does not restore it. 1CRN uses the same frame and its unenforced N–CD edge happens
to lie inside the interval. This is not evidence of a closure-aware constructor.

## Why the larger masks do not contradict that diagnosis

In all six CG/CD and CB/CG/CD cases, availability and deterministic first-valid
graph ordering select **C–CA–N -> CD**, followed by **C–CA–CB -> CG** (with the
existing stereo cone for CB if CB is also generated). N–CD is now an explicitly
constructed length edge, but **CG–CD** is the generated-involving ring edge which
is not enforced as a placement length. Those six saved original proposals pass
their old audit, not because every ring edge was jointly solved, but because the
remaining edge and all assessed terms happen to be within bounds for those poses.
No frame was changed or retried in this diagnosis.

This identifies a **mask-dependent tree-frame limitation on a cyclic graph**.
It is not loss of the bond from chemical authority, and not a reason to remove or
weaken the validator's ring-edge check.

## Fixed observations: what is and is not established

All three CD-only cases pass the conservative necessary distance-shell and
cosine-interval checks for the two observed neighbours CG/N and the authority
CG–CD/N–CD lengths plus angle at CD. The observed CG–N distances are:

- 1TIM: **2.397542283256 A**;
- 2CI2: **2.327282965176 A**;
- 1CRN: **2.388517113190 A**.

These outward-rounded scalar tests do NOT prove joint feasibility of every angle,
plane, stereo and context constraint. No sphere-intersection points, alternative
poses, optimizations or new candidates were constructed.

There is also stronger, explicitly **post-hoc** evidence for the OLD own-request
scope: the retained native-reference additions pass its complete independent
heavy audit for all nine cells, including both rejected CD-only cells, with every
known observation unchanged. They are witnesses for that assessed original scope,
NOT constructor inputs, new proposals, full-complex chemistry or native-H proof.
Thus the two failed chosen poses do not demonstrate that the fixed observations
must move. No native coordinates were used to choose a new frame or torsion.

## Original-own versus COMMON v2-port geometry

The original-own rejection must not be conflated with the separately declared
v2 real incoming amide-port constraints. The original proposal never received
those added constraints during construction.

For the raw CD-only proposals, COMMON geometry adds these plane deviations from
the observed N/incoming-C/CA plane, versus the unchanged **0.05 A** bound:

| Source | Original CD plane distance (A) | Retained native CD plane distance (A) |
|---|---:|---:|
| 1TIM | 0.296695337007 | 0.197977249771 |
| 2CI2 | 0.807229305594 | 0.032864702967 (passes) |
| 1CRN | 0.776623736968 | 0.237171501348 |

1CRN's original COMMON CD–N–CA angle also measures **98.702217 degrees**, below
the added interval **[106, 130] degrees**. Its original OWN audit still passes.
1TIM's retained native reference fails the COMMON plane bound; 2CI2's reference
passes COMMON. Therefore the native witness for OLD own geometry cannot be
silently promoted to a COMMON-port witness in every system. These are scoped
authority/bound compatibility observations, not a diagnosis that deposited native
structures are chemically invalid. No bound is loosened.

Saved PDBFixer coordinates are diagnostic comparators, not a universal escape:
1TIM-CD violates its own target bond-length bounds at both CG–CD and N–CD; 2CI2-CD
violates the own CA–N–CD angle bound. All complete saved arms/references are
measured under both scopes; no failed external candidate becomes a projection
parent or a new accepted structure.

## Verification and retained evidence

Private report: `tmp/chemical-reconstruction-closure-diagnosis-v1/analysis-01/diagnosis.json`.
It binds all nine saved evidence hashes and original/common request, original plan
and raw/rounded candidate digests. It contains every assessed target bond/angle/
explicit-constraint row, fixed/generated roles, frame targets/residuals, omitted
length edges, all available comparator/reference tables and necessary-support
intervals. There are no unexplained or scalar-only audit findings in the actual
nine-cell diagnostics. Unit tests keep admission-short-circuit scalar observations
separate from the independent audit's emitted findings.

Both execution and saved-report reproduction disable 29 known coordinate entry
points, including graph constructors and `extend`/`two_angle_tip`, PDBFixer/shared
rebuild, projection/H builders, optimizers/search and rotations. The existing
WP2P source/authority/template/global-L/H/native readback runs first. Metadata-only
Topology is permitted for its existing global-L detector. **Zero new coordinate
proposals**, zero changed old inputs/results; the old 9/7/2 denominator is unchanged.
Known-entry-point guards are not a general Python call-graph proof.

See the [validation record](../plans/chemical-reconstruction-closure-diagnosis-validation.json)
for exact artifacts, hashes and test/typing checks. Prior WP2O native worsenings,
all v1/altloc refusals, WP2F timeout and broader scientific/production gates remain.

## Next decision, not an implemented fix

A subsequent prospective slice should declare a **generic cycle-aware generated-
subgraph construction** contract: jointly enforce every generated-involving ring
edge and the relevant angular/stereo/port/context constraints while leaving all
observations fixed. It must support full missing subsets, not special-case PRO/CD
or substitute a source/template/native pose. It needs its own raw/rounded audits,
refusal evidence and independent mask/complex validation.

Do not silently let these rejected baselines through the existing projector gate,
relabel them as passing parents, switch to another frame after seeing native RMSD,
or retroactively replace WP2P outputs. Choosing such a new construction/admission
policy requires a separately frozen protocol before new coordinates are proposed.
