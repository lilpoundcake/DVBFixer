# cc89 local refinement diffusion experiment

Status: development experiment requested 2026-10-23; no production promotion.
Input: supplied 9CQD derivative, frozen seed-7 request and immutable source.
See [PATCHR comparison](../research/patchr-gap-reconstruction-comparison.md).
Completed results: [9CQD LRD ablation](../research/cc89-lrd-ablation-9cqd.md).
The original ten candidates fail independent acceptance with geometry defects.
The six generated+anchor minimization additions pass the tested geometry checks,
but still fail strict original-coordinate preservation. No production promotion.

## Experiment matrix

- Initial generation: native conditioning versus per-step fixed-atom reinjection.
- LRD mobility: generated-side boundary atoms only versus both generated and
  observed atoms within a ±3-target-residue boundary window.
- Processing: OpenMM-only, LRD-only, LRD followed by OpenMM.
- Ten distinct final outputs: two OpenMM controls and eight LRD-containing
  outputs. The OpenMM-only row is shared between the two LRD mobility settings.
- Regenerate both initial arms locally with identical model/runtime and seed;
  do not use the historical archive as the only non-reinjection control.

"Reinjection" here is the initial 500-step sampling factor, as in the preceding
experiment. LRD always resets atoms outside its active window after each step:
this defines locality, rather than constituting an additional experimental arm.
LRD is a separate 25-step low-noise cc89 denoising loop (sigma 1.5 Å to 0,
rho 7, positive minimum 0.0064 Å, Euler step scale 1.2, no extra churn).
This is PATCHR-inspired inference with cc89, not a numerical reproduction of
PATCHR's Boltz/Protenix implementation. Full candidate sequence is fixed; no MPNN.
Noise seed is fixed before evaluation and shared at common atoms across mobility
arms. All original context/physical-position/neutral-phase policies are retained.
Both mobility arms use the same centroid of original represented observations;
the frame must not depend on the chosen mobility mask. A preliminary run with
mask-dependent centering is retained as a confounded pilot, not the final paired
comparison. Its frozen initial candidates are reused unchanged for the corrected
frame comparison.

## Acceptance and evidence

- [x] Unit-test boundary identities, mask semantics, schedule and hard projection.
- [x] Freeze inputs, engine/adapter/checkpoint/config hashes and settings.
- [x] Generate both initial arms with complete per-gap callback evidence.
- [x] Run both LRD mobility arms on each initial candidate.
- [x] Run identical CPU OpenMM v8 controls/combinations, without solvent.
- [x] Independently evaluate every candidate against the ORIGINAL source/request.
- [x] Additionally inspect the whole mobile window and its outer junctions;
  the normal generated-only geometry neighborhood is insufficient when source
  residues are intentionally released.
- [x] Report strict acceptance, local geometry and observed-coordinate displacement
  separately. Do not restore mobile source coordinates silently to make a pass.
- [x] Record results, failed attempts, outputs and qualification limits.

Observed-mobile LRD intentionally tests a different scientific constraint.
The normal validator must still reject excessive source displacement. Geometry
metrics may explain an improvement but cannot turn such a result into accepted
production output. OpenMM following LRD continues to move generated atoms only;
it freezes observed atoms at that LRD candidate's positions, not at hidden
restored original positions. A failed refinement is retained as failure evidence.

## Extension: generated atoms plus immediate anchors

User-requested additional factor: repeat each OpenMM-containing branch with
the immediate left/right anchor residue heavy atoms mobile, alongside generated
atoms. For 9CQD these anchors are W25/W29 and X139/X142, not the whole ±3 LRD
window. Keep seed 7, CPU v8, 500 iterations, eight restarts and 0.25 Å
perturbation unchanged; keep the original generation and LRD outputs unchanged.

- [x] Add/test the explicit anchor identity mask, including insertion codes.
- [x] Run six additional minimizations (two initial arms × three input types).
- [x] Independently check original-source drift, anchor/non-anchor drift, full
  mobile neighborhoods and their outer junctions.
- [x] Keep old comparison/PDBs intact and add new labelled visual-review copies.
- [x] Record scientific outcomes without loosening production acceptance.

`scripts/run_cc89_lrd_ablation.py evaluate --minimize-anchors` writes
`processed-with-anchors/` and `comparison-with-anchors.json` separately.
The original request's generated atom/residue mask is never redefined for
validation: the enlarged mask is passed only to the scientific minimizer.
Observed-mobile LRD residues outside the immediate anchors remain frozen at
their LRD positions during this extension. Releasing the whole original window
in OpenMM is not implied by this request.

## Deferred: joint reconstruction of nearby gaps

Status: **deferred at user request — no suitable test structure yet**.
Do not implement or benchmark joint gap generation in this experiment.

When a suitable structure is available:

1. Identify interacting missing regions from observed source-only evidence;
   do not inspect a reference to choose favorable gaps.
2. Freeze sequence, masks, context, seeds and independent versus joint arms.
3. Compare joint sampling against immutable-source independent sampling with the
   same coordinate-preservation, geometry and cross-gap clash checks.
4. Require an explicit policy for any overlap between mobile boundary windows.
5. Keep accuracy, geometry and source motion separate; no reference-guided seed
   selection or general promotion from a single favorable example.
