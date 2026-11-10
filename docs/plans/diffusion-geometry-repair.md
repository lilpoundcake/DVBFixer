# Diffusion geometry repair: implementation and follow-up

Date: 2026-11-10. P0/P1 implemented; P2 private workflow and anchor-calibration
pilots executed after approval. P2 product integration and P3 remain open.
MODELLER remains the default backend. No additional neural sampling is approved
or needed for the initial retained-coordinate controls.

## Evidence and invariants

The ten-complex, five-seed cc94 MPS cohort produced 27 accepted and 23 rejected
candidates. Nine cases used the original code; 1N8Z was a separately versioned
input-preparation fix follow-up. This is neither an accuracy benchmark nor broad
qualification. Retained evidence is below the original checkout's
`tmp/real-complex-cc94-2026-11-06/`, including the dated failure analysis.

Twenty-one rejected candidates retained geometry defects; two 1JPS candidates
had clean geometry but exceeded the observed-anchor movement bound. All 50 raw
candidates required repair. Native/PDB comparison covered 275 gap samples:
internal distances were preserved to PDB rounding, while source restoration
moved 59 initially in-range junctions outside the 1.20–1.45 Å interval. These are
raw junction findings, not 59 final rejected candidates.

Keep original inputs, raw proposals and original outcomes immutable. Every new
repair trial is separate post-hoc evidence, never a replacement cohort result.
Preserve generated masks, atom identities, template ownership, disulfides,
unselected physical records, and zero D-Cα. Never silently rewrite frozen requests.
New numerical algorithms need an explicitly request-bound version; the final
independent validator and current acceptance limits remain unchanged.

## P0 — match search and acceptance (approved)

- [x] Use the full independent geometry assessment for candidate ranking and
  early exit, including Ramachandran, χ1/χ2 and explicit junction connectivity.
- [x] Include per-anchor movement eligibility in state selection and add bounded
  positional protection without relaxing the final limits.
- [x] Preserve historical request decoding and algorithm selection; retain
  rejected trial states and useful search diagnostics.
- [x] Test 1GC1/101, 9CQD/102 and 1JPS/101/103 from retained coordinates, alongside
  previously accepted controls. Do not tune thresholds to these cases.

## P1 — geometry before steric relaxation (approved)

- [x] Add generated/mobile-only projection toward force-field equilibrium bond,
  angle and ring geometry, with source atoms fixed and chirality protected.
- [x] Optimize generated junctions against the fixed observations before gradually
  introducing steric interactions; finish under the full force field.
- [x] Use finite bounded stages, retain diagnostics, and independently validate
  all outputs. Do not freeze the distorted raw bond lengths as constraints.
- [x] Exercise compressed bonds, broken rings, junction defects, exact fixed
  observations, disulfides, chain case and insertion codes in regressions.
- [x] Run unchanged retained-coordinate controls without claiming new sampling,
  structural accuracy or broad qualification.

### Implementation and evidence status

New requests bind `geometry-first-bounded-anchors-v2`; decoded
`geometry-gated-anchors-v1` requests keep the old numerical route and limits.
No-policy requests also retain their historical route. V2 ranks the actual
PDB-rounded projection by admission failures, per-anchor bound violations, then
all actionable full-validator findings; full-force-field energy breaks ties.
An all-zero score is required for successful early exit. This is search ranking,
not a replacement for independent publication validation.

`geometry_conditioning.py` introduces finite bonded-only (up to 100 iterations),
two soft-repulsion (up to 75 each), and full-force-field phases. Force-field
equilibrium bonded terms optimize rings and junctions; this is Cartesian
conditioning, not a guaranteed exact inverse-kinematics closure solver. The
original full nonbonded interactions are restored for the final phase. Temporary
FF hydrogens are discarded on projection. One-sided heavy-atom/Cα and grouped
backbone RMSD walls protect existing movement bounds; the independent audit
still rejects any overshoot. Shared 1:1 SG pairing is applied before snapshotting
temporary v2 disulfides, so an extra distance-inferred bond cannot be promoted
to a mandatory second SG partner.

Private trials retain PDB projections, hashes, full assessments and movement
reports; low-level diagnostics identify phases and scores. Original inputs,
raw proposals and cohort outcomes remain unchanged. Replays are separately
labelled post-hoc repair, not new neural samples.

Initial Reference tests cover equilibrium bond projection, compressed proline
ring closure, fixed flanks, full-score-driven restart selection, force-group
switching and anchor-wall gradients with case-sensitive/insertion-code identity.
The seven declared retained-coordinate controls completed on CPU, without new
sampling. Independent full validation passed six of seven:

| Retained proposal | Original outcome (unchanged) | V2 repair validation |
|---|---|---|
| 1GC1 / 101 | rejected | passed, generated-only |
| 9CQD / 102 | rejected | passed, generated-and-anchors |
| 1JPS / 101 | rejected for anchor drift | passed, generated-only |
| 1JPS / 103 | rejected for anchor drift | passed, generated-only |
| 1KTZ / 7 | rejected | still rejected |
| 1BRS / 7 | accepted | passed, generated-only |
| 1G9M / 7 | accepted | passed, generated-only |

Both 1JPS controls now pass without moving any observed anchors. In 9CQD the two
selected anchors W25/W29 have backbone RMSD 0.1769/0.4749 Å, Cα displacement
0.0768/0.2918 Å and maximum heavy displacement 0.3013/1.0731 Å. No ceiling was
raised. These post-hoc controls are not a new cohort acceptance rate or evidence
of structural accuracy. The original 27/50 outcomes remain unchanged.

1KTZ/7 retains 21 severe clash findings, concentrated around generated A6 and
C4/C5, plus C4–C5 C–N = 2.1689 Å and a bad C4 angle. Bonded preconditioning does
not guarantee that a deeply interpenetrating loop escapes its starting basin;
the final full-force-field phase can still trade bond strain against sterics.
This negative result is retained and rejected, not labelled repaired. Exact
closure and a more robust conformational search remain additional work; P2/P3
are still plans, not hidden fallback capabilities.

Evidence and exact software snapshots are under
`tmp/geometry-repair-v2/controls-01/`: `declaration.json`, per-case
`protocol.json`, unchanged raw artifacts, trial projections, phase diagnostics,
full validation and independent custody/readback records. All six accepted
public bundles passed independent readback. Prepared-source
publication retains the original repair bytes/report and frozen snapshot, then
explicitly binds a new request snapshot containing only the post-hoc refinement
policy change; it does not rewrite historical preparation evidence. The final
`independent-audit.json` SHA-256 is
`f6231f7ac9cc29afc8c1023de244d18f1786d99c3850b3f9ae35a44abfe42395`.
The full fast suite
passed 2,222 tests (three existing optional Torch/Biotite skips; 33 slow tests
deselected); the final focused set passed 48, including one later-added soft-core
force test. Ruff and agent-document maps passed. Slow integrations and GUI
testing were not rerun for this Python-only numerical change.

## P2 — per-gap portfolios and repair scheduling (private implementation)

2026-11-10 continuation: a separate private multi-parent contract, bounded beam,
issue-local generated-only refinement and independent lineage/mobility audit
are implemented. A retained 1KTZ portfolio plus D131–136-only repair passed full
geometry/custody; the original 0/5 remains unchanged. Public bundle/`.dat` and
CLI/Homology/GUI integration remain open. See the
[evidence and remaining gates](../research/gap-portfolios-and-anchor-calibration.md).

The 1KTZ retained stage-1 reports contain at least one locally clean option for
each of its eight gaps, but not in one shared-seed candidate. This is potential
coverage, not proof that a recombined complex will pass.

- [x] Define per-gap diagnostics and issue-local repair budgets, preserving good
  generated regions instead of perturbing every gap on every retry.
- [x] Design explicit multi-parent lineage for gap choices from different seeds;
  do not mislabel a recombination as a single original seed or mutate schema-4
  frozen requests/public bundles to hide the additional selection operation.
- [x] Select and materialize only in Python; no scientific recombination in Node
  or the GUI. Preserve common source/frame and template-ownership bindings.
- [x] Reassess cross-gap clashes, junctions, global chirality, metadata and final
  whole-complex acceptance after assembly. Local cleanliness is not publication
  authority, and original rejected candidates remain rejected.
- [ ] Predeclare a separate evaluation and resource budget before new sampling.

## P3 — explicit alternative reconstruction (planned, not implemented)

1HE8 A322–352 remains problematic in all five candidates. Retrying one Cartesian
proposal indefinitely is not an adequate general reconstruction strategy.

- [ ] Specify an explicit opt-in alternative/loop-search policy for failed gaps,
  potentially using existing MODELLER infrastructure. No silent fallback when
  the user selected diffusion, no invented successful sampling provenance.
- [ ] Define eligible chemistry, fixed-observation or bounded-anchor contracts,
  fragment/complex handling, licensing/readiness checks and failure retention.
- [ ] Bind each alternative proposal and subsequent refinement to its true
  source, method, parameters and generated atom ownership.
- [ ] Require the same final geometry, identity, chirality and publication gates.
- [ ] Evaluate independently of the original cc94 cohort and retain negative
  outcomes; do not replace the baseline acceptance rate with post-hoc rescues.

## Anchor-limit calibration (pilots executed, limits unchanged)

Predeclared ubiquitin/NMR and lysozyme comparisons have executed. Leave-anchor-out
local/global frames, exact identity, exclusions, contacts and source QC are
reported separately. Initial 1LYZ has source defects and is retained as a QC
warning; a separately declared 1AKI/193L control followed. These two-protein
methodological pilots do not establish a production ceiling. See the
[measurements and limitations](../research/gap-portfolios-and-anchor-calibration.md).

The current per-anchor limits—backbone RMSD 0.5 Å, Cα displacement 0.5 Å, maximum
heavy-atom displacement 1.25 Å—are engineering bounds, not measurements of
experimental uncertainty. The user requested investigation of a higher ceiling;
this is not permission to raise it merely to accept existing failures.

- [x] Predeclare pilot candidate bounds and cases independently of comparison outcomes.
- [ ] Compare independent experimental models, alternate conformations and
  repeated/refined structures in a common local reference frame. Account for
  experimental method/resolution, occupancy and model uncertainty.
- [ ] Treat B factors as supporting evidence, not a direct displacement radius
  or an automatically calibrated error estimate.
- [ ] Measure changes in backbone, sidechains, interfaces, covalent links,
  density support where available and unselected observations; distinguish
  relative domain motion from local coordinate correction.
- [ ] Report sensitivity at existing and proposed bounds, including geometrically
  clean candidates currently rejected only for drift and damage introduced by
  excessive mobility. Keep case discovery and prospective validation separate.
- [ ] Any adopted increase requires a new request-bound policy version, explicit
  movement reports and independent revalidation. Historical requests retain
  their original limits and serialization.

## Tracking

See the [canonical checklist](implementation-checklist.md),
[anchor-relaxation handoff](geometry-gated-anchor-relaxation-handoff.md), and
[plan index](README.md). Private P2 and calibration research are not public shipped
capabilities; P3 remains planned only.
