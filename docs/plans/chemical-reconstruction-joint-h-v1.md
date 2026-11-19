# WP2G: joint generated-only sp3 H-group orientation

Date: 2026-11-19. Prospective protocol, baseline `8daa3b6`; persistent
`feature/chemical-reconstruction`. Freeze before any joint-fit numerical controls.

## Bounded change and documentation plan

Retain the complete WP2F H inventory, source preparation/audit, chemical states,
canonical names, all observed coordinates/occupancies/context, constraints and
independent geometry/CIP validators. Add a separate digest-tagged private plan
for joint neutral-carbon sp3 groups only: three observed heavy neighbours/one H
or two observed heavy neighbours/two H. Replace their independent cone proposals
with one shared proper rotation per group. Do not change backbone amide,
aromatic CH, hydroxyl or terminal-three-H constructors/phases.

Documentation: this protocol, joint-H validation JSON, active handoff and agent
task map. Pending links: none. Old protocols, results and old default plan digests
remain unchanged. No public/default/dependency/data-license capability follows.

## Objective and fixed algorithm

For each group, subtract its fixed parent centre and normalize every source ideal
and observed heavy-neighbour vector separately. Minimize the equally weighted sum
of squared directional chord residuals with a SINGLE determinant-+1 rotation.
No translation, scale, observed-atom movement, group-specific fitted weights,
phase search, H relabelling, stochastic seed or iterative optimizer is allowed.

Reuse the existing `model.diffusion.geometry.weighted_kabsch` proper-rotation
owner on antipodal clouds `[v1,-v1,v2,-v2,...]` and matching target clouds. Their
centroids are mathematically zero and the objective is twice the desired
origin-constrained directional objective. This is reuse of a geometry helper,
NOT a diffusion request, sampler, source repair or frozen-profile intervention.
Ignore/reject non-negligible numerical translation rather than applying it.
Apply the shared rotation only to source ideal H directions, use each pinned
source bond length, and re-anchor at the UNCHANGED observed parent coordinate.
Thus H-H angular geometry is transported jointly, while incompatible heavy/H
geometry can still be rejected by the unchanged independent validator.

Require finite nonzero directions and a non-collinear pair among ALL heavy
neighbours in both frames: maximum cross-product norm must exceed 1e-8. Require
proper orthogonal output (`R.T R = I`, determinant +1 within absolute 1e-10)
and translation norm below 1e-10. No alternative pair/frame is selected on failure.
All parent/heavy anchors must be unit-occupancy observations, every generated H
belongs to exactly one complete admitted group and has a single bond to its
neutral-carbon parent. No external/covalent/H ownership guessing.

Represent group membership/all heavy identities explicitly in the new plan;
keep per-H placement coverage exact and all original ordered-H/full-ligand
stereo constraints. Constructor caches one fit for the whole group. Independent
candidate audit validates group authority/coverage and geometry without calling
the fitter; constructor success never establishes acceptance. Existing full-H
conditional CIP and its H-coordinate/pro-R/pro-S limitations remain explicit.

## Immutable validation policy

UNCHANGED WP2F bounds: observed-heavy 0.10 A / 12 degrees; generated-involving
bonds/angles 0.02 A / 5 degrees; H planes 0.005 A, heavy planes 0.05 A; minimum
separation 0.7 A; normalized stereo minimum 0.05. Keep all fixed-heavy and complete
four-ligand checks. Round only additions, independently audit raw/rounded and
retain every construction, geometry, source or CIP refusal. No acceptance based
solely on directional-fit RMS residual, graph valence or CIP label.

## Twelve declared comparison entries

Reuse EXACT WP2F table, full masks and source/authority pins, as an explicitly
retrospective case comparison with a prospectively frozen new algorithm:

- Five nominal synthetic ALA/THR/ILE/LYS/PHE controls at incoming angle 121.96;
  one synthetic ALA incoming angle 100 control.
- Four targets in 1UBQ: first interior ALA28, ILE3 and LYS6 by file order and
  exact PHE45; first interior 1FR2 THR11; first interior 1A3N ALA.

Twelve entries, not twelve independent folds. Hypotheses: nominal synthetic and
prior ALA/THR observed successes remain; joint orientation may resolve ILE/LYS
angular deficits, but NO required pass count is asserted. PHE fixed-heavy ring
and 1A3N altloc refusals remain; distorted amide remains rejected because its
strategy is unchanged. Preserve the original WP2F ILE/LYS failures irrespective
of new results. Do not broaden masks, alter pH/state, choose another target or
widen bounds after execution. No native-H reference scoring or conformer ranking.

Check frozen protocol and CCD/FF/naming/observed pins before exclusive output
reservation. Retain original/prepared bytes, old view, new plan/group ledger and
raw/rounded results plus conditional CIP. Independently read back saved files
and remeasure without proposal construction. Replay old control artifacts byte-
for-byte as regressions only; do not count repeats as scientific evidence.

## Regression and remaining gates

Test proper SO(3) and rigid equivariance, rank/degeneracy refusals, all-heavy use,
group H-H/length preservation, exact membership/coverage, missing or unknown
group plans, wrong element/charge/bonds, observed/generated ownership, incomplete
groups, unchanged old defaults, corrupt candidate and context-clash rejection.
Run focused fast tests, separate seeded/renumber boundaries, Ruff, scoped strict
mypy and full fast suite. Preserve WP2F's unresolved large combined slow-run
timeouts; no skip/assertion/platform changes or false combined-pass claim.

This remains case-specific private complete-selected-H research. Measured/partial
H, all20/linked chemistry, stereospecific H-name qualification, actual-neighbour/
polymer chemistry, native-H/hydrogen-bond quality, broad complex qualification and
production/public/default/dependency/data-license gates remain OPEN.
