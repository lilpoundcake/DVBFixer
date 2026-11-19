# WP2D: explicit metric geometry and external-anchor amide H

Date: 2026-11-19. Prospective protocol, baseline `1ae0de9`; continue only on
`feature/chemical-reconstruction`. Freeze before new metric controls execute.

## Plan and authority

Implement an extensible private metric-geometry plan alongside the unchanged
template-only plan. Explicit, provenance-bound bond/angle targets may replace
template measurements only for geometry involving generated atoms. Do not
silently alter CCD coordinates/chemistry or reinterpret free-amino-acid N-H
geometry as amide geometry. Existing plan digests/default behavior must remain
unchanged when no metric plan is used. Constructor and measurement-only auditor
must independently consume the same bound authority; the auditor never calls
the constructor. Retain raw and generated-only PDB-rounded results, including
failures. No optimizer, alternative-engine retry or observation movement.

Reuse shared requests/ownership, conditional peptide/naming proof, existing
geometry math and audit/clash infrastructure. The first metric-driven placement
is planar three-coordinate geometry with a local centre, two observed bonded
neighbours (one may be an explicit covalent external connection) and one generated
atom. Three ordered anchors and all required length/angle metrics are mandatory.
This is an extensible geometry boundary, not a single-H product or a replacement
for complete heavy-sidechain/H objectives. General cycles and tetrahedral H
groups remain separate strategies.

Before this freeze, planning inspected only the existing ff19SB H-N bond and
C-N-H / XC-N-H rows: 0.101 nm, 120 degrees and 118.04 degrees, with equal angle
stiffness 418.4. No new coordinates or numerical controls were executed. Pin
`FF/amber19sb.ff/ffbonded.itp` SHA-256
`b6fa511c5bd3f11a19484a418d10e32972e2c19788d464466926c7d8fbde1034`.
Reuse the existing RTP parser/inventory pins and canonical naming owner. A narrow
parameter adapter selects unique exact function/type rows, retains exact tokens,
source pin and declared units; no wildcard resolution or partial-charge-derived
formal charge. The RTP must explicitly identify target N/H/CA as N/H/XC. Refuse
wrong types, duplicate/absent rows, unequal stiffness or changed pins. Force-field
equilibrium metrics are an explicit research geometry authority, not a calibrated
experimental geometry library or parameterization of the retained complex.

## Construction and independent audit

For observed neighbour angle alpha, planar exterior H angles must sum to
360-alpha. Use the analytical equal-stiffness solution:
theta_first = (360-alpha + target_first-target_second)/2. Place the new atom in
the observed neighbour plane on the exterior side, at the authoritative bond
length. This is a declared angular compromise, not exact satisfaction of both
equilibrium angles for every distorted input. Degenerate/collinear anchors,
interior solutions and unavailable chemistry/metric/connection authority refuse.
Measure both angles, length and planarity independently against declared bounds;
do not accept solely because the constructor completed.

Bounds: bond deviation 0.02 A, angle deviation 5 degrees, generated/context
minimum separation 0.7 A, planar tolerance 0.005 A. Observations remain exactly
fixed; PDB rounding applies only to additions. Audit the complete retained request
using existing stereo/planar constraints and graph-distance clash exclusions.
Explicit external links authorize anchor connectivity, not invented inter-template
bond lengths. Never apply missing-template external angle defaults.

## Six prospective controls

Five positive hypotheses: ALA, THR, ILE, LYS and PHE conditional peptide states,
each with only backbone H absent. Other projected ideal CCD atoms remain observed.
The incoming C is an explicitly SYNTHETIC opaque neighbour positioned from N/CA
and the N/CA/C plane, at 1.335 A and Cprev-N-CA angle 121.96 degrees. No withheld
H coordinate participates in anchor placement or construction. Use one fixed
rigid transform for a second positive orientation in regression tests, not an
independent scientific cohort. Source bytes are a digest-bound synthetic record,
not original PDB observations. Conditional ALA witnesses still do not prove the
actual opaque neighbour chemistry.

Sixth prospective control: ALA with incoming angle 100 degrees, unchanged
authority/bounds; expect retained raw and rounded angle rejection, not widening,
anchor movement, geometry repair or a different parameter choice.

Retain six full request/plan/transport/authority/result records and a summary.
Replay/regression tests are not independent evidence. Negative regressions cover
metric/request binding, missing/duplicate metrics, non-generated overrides,
missing/wrong external links, collinear anchors, bad parameters/pins, a corrupted
H coordinate, changed length/angle, a retained-context clash, old digest stability
and unchanged WP2B H refusal. No WP2B observed source is repaired or rerun.

## Verification and documentation

Run focused heavy/H/geometry and legacy-boundary tests, Ruff, scoped strict mypy,
then full fast suite. Record unchanged transitive legacy typing failures honestly.
Verify all original fixture/FF bytes, old controls/default digests, final source
and report hashes and local links. Documentation plan: this protocol, amide
validation JSON, active handoff and agent task map; pending links: none.

This is private synthetic metric/geometry evidence only: no observed-complex H
qualification, whole-polymer state/CIP, arbitrary partial H/subset coverage,
public/default/dependency replacement or production/data redistribution approval.
Source admission/preparation is still a distinct prospectively specified next
slice, preserving original files, full context and prior refusals.
