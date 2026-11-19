# WP2H: explicit planar-N and thiol H groups

Date: 2026-11-20. Prospective protocol; baseline `357c821`, persistent
`feature/chemical-reconstruction`. Freeze before source acquisition or numerical
controls. This is private research, not all20 qualification or a public capability.

## Scope and documentation plan

Retain old complete-H and joint-H entry points/results unchanged. Introduce a
separate opt-in extended compiler using the same request, naming, authority,
constraints, retained-complex source adapter, metric constructors and validators.
Extend H groups, NOT heavy-atom or source repair, to authority-proven trigonal
nitrogen and neutral thiol sulfur. Keep heavy reconstruction as an equal track;
this slice does not qualify missing heavy sidechains for these new residues.

Documentation: this protocol, independent validation JSON, source fixture
provenance, active handoff and task map. Pending links: none. User's continued
implementation request authorizes this bounded slice; no public/default changes.

## Sources and exact declarations

Download unmodified RCSB CCD CIFs for ASN, GLN, ARG, CYS, TRP, GLY, PRO, HIS from
`https://files.rcsb.org/ligands/download/{ID}.cif` into a separate fixture directory
`tests/fixtures/chemical_reconstruction/ccd_h_groups_v1/`. Retain CC0/source URLs,
bytes, SHA-256, revision and atom/bond counts in its own manifest. Existing seven-
component CCD manifest and all old pins remain unchanged. Pin and commit the new
source manifest before numerical construction controls. No PyRosetta data.

Use the existing neutral-ALA conditional chemical witness and exact AMBER19SB
STANDARD RTP inventory for the SAME named state. HIS must not silently become
HID/HIE/HIP; PRO must not obtain a new leaving-atom projection by inference.
Missing inventory, failed state proof or source/geometry admission is a retained
refusal, not permission to rename, protonate, remove atoms or choose another case.

Sixteen entries: each of the eight named states above once synthetic and once
observed. Synthetic heavy coordinates are the pinned local ideal target plus a
real-in-request incoming carbon arranged at Cincoming-N-CA angle 121.96 degrees,
N-Cincoming distance 1.335 A, with no other generated atoms or real neighbours
replaced by witnesses. The ALA chemical witness is virtual chemistry only.
Observed declarations: first file-order interior matching residue in chain A of
1UBQ for ASN, GLN, ARG, GLY, PRO and HIS; chain A of 1FR2 for CYS and TRP. Preserve
all source records except the already explicit dangling-CONECT preparation policy,
all observations, links and the unchanged opaque fractional-context policy.
No target substitution, altloc selection, occupancy rewrite or hidden source fix.

## Fixed construction strategy

An N with three local graph neighbours is admitted only if the exact parent-plus-
neighbour set occurs in the existing pinned CCD authority's SP2 planarity groups.
This reuses the existing sanitized explicit graph and does not infer protonation
or hybridization from residue name or coordinates. Require supported N charge
0 or +1, H single bonds, and exact one-H/two-heavy or two-H/one-heavy inventory.

- One-H/two-heavy planar N: use the existing metric planar exterior constructor,
  source N-H length and both ideal heavy-N-H angles, plus the unchanged strict
  generated-H plane constraint. Both heavy directions participate.
- Two-H/one-heavy planar N: use existing chain internal coordinates with one
  deterministic heavy anchor adjacent to the heavy neighbour, selected by exact
  lexical order among nondegenerate ideal frames. Both H retain their pinned
  named torsions and lengths. Validate H-H angle and N/three-neighbour planarity;
  also require the heavy anchor and complete N-H2 group plane at 0.005 A.
- Neutral, nonaromatic S with exactly one single-bond heavy neighbour and one
  single-bond H: reuse the pinned heavy-chain torsion strategy, retaining its
  named ideal phase. This is thiol only, not disulfide/thiolate chemistry.

All existing groups and joint neutral-carbon fitting are reused unchanged. No
phase search, torsion optimization, H relabelling, new geometry engine, observed
motion, chemistry edit or fallback. Mark the extended origin policy separately.
Old compiler still refuses previously unsupported groups.

## Immutable evaluation and regression policy

Same bounds: fixed-heavy 0.10 A / 12 degrees; generated 0.02 A / 5 degrees;
separation 0.7 A; H planes 0.005 A; heavy planes 0.05 A; stereo minimum 0.05.
Raw and rounded additions must independently pass the existing geometry and
conditional-CIP audits; no pass count is required. Chemical refusals, invalid
fixed planes, context clashes and constructor/rounded failures remain visible.
Conditional CIP is not H-coordinate or stereospecific-name qualification.

Retain sources, inventory/chemical proof, input request, group ledger, plan and
all results. Independent saved readback remeasures without generating proposals.
Test exact graph/SP2 authority admission, nonplanar/tetrahedral-N rejection, charge,
bond/ownership/inventory refusal, named NH2 geometry, SH torsion/length, generated
plane corruption, rigid equivariance and old-default compatibility. Replay old
controls byte-for-byte; do not count replays as independent science. Run focused
tests, split slow boundaries, Ruff, scoped strict typing and full fast tests.

## Limits

Same-authority synthetic consistency is not native-H holdout quality. Observed
cases are bounded retained-complex controls, not independent folds or all20
coverage. Full polymer/actual-neighbour chemistry, measured/partial H, stereospecific
names, orientation/hydrogen-bond quality, rings/heavy repair and broader chemical
states remain open. Preserve WP2F's unresolved combined slow-run timeout. No
production/public/default/dependency or licensed-data gate is closed here.
