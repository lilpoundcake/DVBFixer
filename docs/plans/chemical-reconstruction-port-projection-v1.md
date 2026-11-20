# WP2M: real planar-peptide-port constraints and generated-only metric projection

Date: 2026-11-20; baseline `65fce4d`; persistent
`feature/chemical-reconstruction`. Freeze BEFORE new projection controls.

## Scope and documentation plan

A separate private generated-heavy proposal, not a retry inside WP2L. Reuse
canonical identities/graphs/requests, source transport, all existing heavy audit
and constraint measurements. Add a generic bounded generated-coordinate metric
projector plus a narrow explicitly typed tertiary peptide-N port declaration.
Do not implement PRO coordinates with a per-residue positional patch. Full and
partial generated subsets are supported by the generic core, subject to bounds;
unqualified chemistry remains refused, never removed from retained context.

Documentation: this protocol, validation JSON, active handoff and task discovery;
local QA quality note. Pending links: none. User's continued-work instruction
authorizes this bounded same-branch slice. No new empirical/licensed data, public
capability, default/dependency replacement or production permission.

WP2L's PRO19 parent and failed H stage stay immutable historical evidence. A new
port-constrained request has its own digest and explicit source-original lineage.
Neither a passing projection nor a metric fit is complete heavy+H qualification.
Projection-to-H conditioning/publication is a separate integration gate; never
cast a projection into a fake WP2K torsion-grid parent to enter WP2L.

## Geometry authority and actual observed port

Use the unchanged pinned ff19SB `FF/amber19sb.ff/ffbonded.itp` (SHA256
`b6fa511c5bd3f11a19484a418d10e32972e2c19788d464466926c7d8fbde1034`).
Select exact rows only, no wildcard/nearest-type fallback:

- C–N bond: 0.1335 nm (1.335 A).
- C–N–CT and C–N–XC angles: 121.9 degrees.
- CT–N–XC angle: 118 degrees.

The narrow adapter explicitly declares the ff19SB types C/N/CT/XC for real incoming
carbonyl C, target N, target CD and target CA. Verify the target inventory/projection
and exact graph bonds: neutral nonaromatic N, two singly bonded internal heavy
carbon neighbours, no attached H/D; exact fixed unit-occupancy incoming C/N/CA and
declared order-1 covalent incoming C–N link. The CD must be generated. Missing or
wrong typing/declaration fails closed; this is NOT atom-type or protonation inference.
Driver pins RTP and binds the selected target inventory. The real incoming point
is from the original source, never virtual ALA coordinates. The chemical witness
still uses ALA and does NOT establish actual-neighbour/full-polymer chemistry.

Append explicit constraints to the original heavy source request BEFORE projection:
incoming C–N distance 1.335 +/- 0.10 A; the three sourced angles +/-12 degrees;
N/incoming-C/CA/CD plane bound 0.05 A. Keep the original constraints intact as a
prefix, including all ring/stereo and peptide-distance constraints. Do not replace
CCD internal metrics with FF angles: both apply, and incompatibility is a refusal.
Preserve all source bytes/observations, components, state, links, pH and heavy mask.
Prepared-source policy/origin remains genuine, with explicit added port/parent
lineage settings; ordinary request digests and old defaults do not change.

## Generic bounded coordinate proposal

Input: explicit source-bound heavy request, a complete passing rounded baseline,
and an exact ordered noncollinear fixed-anchor frame. No native reference, heldout
coordinates, H outcome, empirical rotamer data or success-dependent settings.
The original baseline must pass its original raw/rounded heavy geometry. Do not
project an unqualified baseline to hide its original failure. New port constraints
can fail on the baseline: that is the declared reason for the NEW proposal.

Only generated heavy coordinates are optimization variables. No observed/mobile/
temporary coordinates, atom addition outside mask, H creation or reflection.
Bound: at most 16 generated atoms/48 variables, one target, complete known heavy
inventory, at most 512 instantiated target metric/constraint residuals; refuse
over-budget inputs instead of truncating or cropping the context. Keep full retained
observations for final geometry/clash auditing. Other constraints must be among
distance/ring-closure/exclusion, angle, plane and ordered stereo (no torsion yet).

One deterministic local solve from the same rounded baseline. Variables are
displacements in the proper observed-anchor frame, not absolute laboratory axes.
Use SciPy least_squares, method `trf`, 2-point Jacobian, linear loss, x_scale=1,
ftol/xtol/gtol=1e-10, max_nfev=500; no random/multistart/phase search or H-dependent
restart. SciPy is optional private research machinery: absence is a refusal, not
permission to alter installation dependencies or silently substitute an algorithm.
Coordinate bounds are +/-2 A per local displacement component; an independent
per-atom Euclidean movement cap of 2 A applies to raw AND rounded proposals.

Residuals: ALL instantiated target heavy bond lengths and local heavy angles from
the unchanged graph authority, normalized by existing 0.10 A /12-degree heavy
tolerances; explicit distance/angle constraints use interval midpoint/halfwidth;
plane signed distances use their declared bound; ordered stereo contributes a
one-sided minimum-volume violation; exclusion uses only a lower-bound violation.
No native fit, positional pull to the heldout conformation, steric/FF energy,
electrostatic/H-bond/solvation term or rank-dependent adjustment. Fixed-only residuals
remain present and uncorrectable. Accept using independent existing raw/rounded
geometry/constraints/clash/stereo audits, movement/ownership/rounding binding AND
solver convergence. Optimizer residual/cost decrease is not acceptance or native
quality. Retain raw/rounded failures and solver disposition even when no pass.

## Frozen controls, readback and limits

Synthetic graph/port controls and corruption tests use tracked authorities and
synthetic observations, not ignored prior outputs. Test generated-only ownership,
proper rigid-frame equivariance, ring closure, complete/partial masks, fixed-source
defects, conflicting constraints, caps/budgets and constructor-free audits. A
distinct projection strategy is rejected by the ordinary single-proposal constructor.
Audit saved coordinates/metrics/constraints/rounding without least_squares or any
coordinate constructor; do not re-solve to validate a proposal.

Prospective source cohort: first file-order interior PRO in chain A of EACH of
`tests/fixtures/protein_context/1UBQ.pdb`,
`tests/fixtures/regressions/1FR2.pdb`, and
`tests/fixtures/protein_context/1A3N.pdb`. Full retained complex and whole sidechain
CB/CG/CD mask, same explicit PRO leaving-H/state authority as WP2L. Selection is
not success-dependent. The 1UBQ case is known-negative continuation; the other
two entries are prospective for THIS method, not new independent protein families
or checkpoint qualification. Refusals remain in the three-entry denominator.
Build a separately labelled unchanged old heavy baseline on each source, then
the new port request and ONE projection; never alter WP2L's saved parent. Extract
exact-name native RMSDs only AFTER projection returns, never as a solver input,
selection criterion or required improvement target. No H-success count is claimed.

Run constructor-free saved readback, byte-identical old controls, focused and split
slow boundaries, Ruff, scoped/source-bound typing and full fast tests. Preserve
ambient legacy typing failure, WP2F timeout and old ASH native-score worsening.
Any later projected-parent H integration needs explicit lineage/admission/audits,
not fabricated observed provenance or relaxation of the 0.05 A port-plane bound.
