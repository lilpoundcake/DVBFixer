# Private hydrogen-initialization guard v1

## Authorization, scope and documentation plan

Continue from `0038b6c` on `feat/hydrogen-initialization-guard`, following the
user's instruction to continue until input is needed. This authorizes a private
prototype and bounded qualification, **not default/public integration, main merge
or release**. The [HB-only result](../research/thr2-hb-only-initializer.md) is local
causal evidence, not qualification of the new selective rule below.

Documentation scope (in order): this prospective contract, a frozen static-cohort
JSON plan, a scalar result note/record, then the canonical handoff update. Track
no future links until their targets exist. Pending links: none. Validate links,
source hashes and absence of coordinate/RNG/state payloads before publication.
The user's autonomous-continuation instruction covers this limited documentation.

## Pure prototype boundary

Keep implementation in `scripts/`, with no import/registration in production.
Accept immutable atom records with `ExactAtomRef` (case-sensitive chain,
insertion code, model, altloc and occurrence), residue name and element, nm
coordinates, undirected bonds and an explicit set of **new H identities**.
Everything outside that set is observed/fixed for the proposal, including input H.
Validate finite bounded triples, identity uniqueness, component residue-name
consistency, unique non-self bonds, endpoint membership and mask membership.
Each new atom must be H with exactly one heavy parent. Unsupported centers are
explicit no-ops, never permission to guess their geometry. No coordinate input
or output is mutated. Refusal is atomic: no partial proposal is returned.

## Prospectively fixed candidate rule

Enroll by graph, not by target name: a canonical-protein carbon (shared canonical
atom inventory and canonical variant map), exactly three heavy neighbors and
one H, all four neighbors in the same component, canonical heavy names, and all
heavy neighbors observed. Only its **new** H may be changed. Other chemistries,
multi-H parents and input H remain unchanged, with explicit accounting.
This is spatial compatibility with the existing heavy tetrahedron, **not CIP**,
L/D correction or a claim that the heavy stereochemistry itself is correct.

For the same deterministic neighbor order, compare
`det(n1-parent,n2-parent,n3-parent)` with `det(n1-H,n2-H,n3-H)`.
Use a dimensionless normalized determinant (divide by the product of the three
vector lengths); magnitude **<=1e-4** is degenerate and refuses for selected
supported centers. Heavy-neighbor distances must be **0.10–0.20 nm inclusive**;
current H radius must be **0.05–0.15 nm inclusive**. Bounds are engineering scope
limits, not force-field equilibrium parameters or general geometry acceptance.

Compatible new H are exact no-ops. For an incompatible H only, use the normalized
sum of **unscaled parent-minus-heavy-neighbor vectors**, preserving its current
parent–H radius. Neighbor-sum length **<=1e-12 nm** refuses. Require the proposed
four-ligand sign to match the unchanged heavy sign, with normalized magnitude
>1e-4. Do not reflect the parent or any heavy atom. The radius-preserving,
selective rule is a new candidate, not the previously tested unconditional
0.1 nm HB intervention. No noise removal for other H, RNG use or template tables.

After all proposals, require every moved H to remain **>=0.05 nm** from every
nonparent atom (including other proposals); otherwise refuse the whole proposal.
This prevents gross overlaps, not all steric clashes, and does not certify the
unmoved initializer. Do not tune bounds after cohort outcomes or retry failures.
Engineering tests cover exact preservation, no-ops, covariance under proper
rotation/translation, permutation invariance, input/observed H, ambiguity,
degeneracy, unsupported chemistry, identity separation and atomic refusal.

## Source adapter and qualification gates

Saved native initialization admission must match atoms/bonds and pre-H observed
identities, validate the H-only System's masses and assert all heavy **and
observed H** pinned, no constraints/virtual sites. Derive the new-H mask from
positive particle masses, cross-check against pre-H absence and element, rather
than trusting a label. Ambiguous older occurrence metadata may be admitted only
when exact keys are unique; never collapse repeated chain occurrences.

Before evaluating the prototype on retained scientific sources, freeze all three
original THR2 initializations plus **all twelve** fixed real-cohort cases/calls,
source hashes and actual route metadata. Static execution only first: no Context,
minimization, new draws or new H sampling. Retain admission refusals, cases without
native calls and unsupported coverage; never relabel requested FF as actual FF.
Original THR2 actual FF XML provenance is missing; keep that state explicit.

Only after static accounting may a separately frozen native experiment be designed
for independent contexts, original physics and stopping settings. Static no-ops
are not dynamic validation. A supported public OpenMM initialization hook has
not been established: no global minimizer monkeypatch may enter production.
The lost original THR2/SER trajectories and earlier cohort failures remain lost/
failed, not replaced by this work. Empirical-table production gates stay closed.
