# Private Mac conditioning continuation: A/C/D

Date: 2026-11-11. User requested continuation until the approved study is finished
or a real user decision is needed. This is private development evidence, never a
public capability or production-default change.

## Frozen next slice

Thirty serial RFD3 calls, six development seeds 7/42/101/102/103/104 per arm:

| Series | New arm | Unchanged evaluation task | Calls |
|---|---|---|---:|
| A | Eight observed sequence flanks | 1UBQ A1–72, hidden A10–13 | 6 |
| A | Flanks plus source-only spatial neighbors within 6 Å of those flanks | Same | 6 |
| C | Backbone-only observed coordinate information, full declared projection | Same | 6 |
| D | Independent gap0 call, excluding the other unknown gap | 1UBQ pair A10–13/A46–49 | 6 |
| D | Independent gap1 call, excluding the other unknown gap | Same pair | 6 |

Reuse existing six-seed full-heavy tiny-gap and jointly sampled pair baselines;
do not resample them or choose seeds from reference scores. This produces 24 new
assembled outcomes: six for each A/C arm and six independent two-donor pairs.
The D comparator changes coupling **and complementary-gap sequence/task size**;
it is not a pure weights-only or fixed-shape intervention. All original observed
residues remain in each D input; the unselected missing region is absent, never
a fixed zero-valued structure or a previous generated donor. Noncontiguous source
segments use explicit `/0` cuts, not invented peptide bonds.

Canonical atom names and complete fixed sequence are transported from the frozen
request. In C, observed sidechains are zero-coordinate unfixed nuisance atoms in
the model, **not newly generated output ownership**. All original observations
are restored exactly in the final source frame. Report represented-observation
native drift and restoration damage separately. Never call this arm strict
source-preserving heavy-atom inpainting.

Before weights load, run the pinned RFD3 input builder on each of the five shapes:
verify every source-identity conditioning mask, fixed sequence and zero unfixed
coordinates. Runtime final-consumer receipts independently check mask counts and
sequence plus six distinct first-noisy arrays per arm. No checkpoint, source,
scientific threshold or old runtime patch is changed.

All 30 calls use the existing explicit MPS seed contract, 200 configured timesteps,
two recycles and one sample. One ML process at a time, pressure-v2, 1800 seconds
per call, RSS8GiB/MPS6GiB, no generic CPU fallback or remote MSA. Maximum declared
sampling wall budget is 15 hours, not a forecast. Existing assets/environments
are reused. Technical failures stop the series and are retained; no implicit
retries. Original 1HE8 outcomes and manual approvals remain unchanged.

After source-only materialization, check the entire original evaluation structure
including observations omitted from local conditioning. Both independently sampled
gaps are minimized together, with exact observations fixed. Only geometry failures
receive the same CPU1thread generated-only 1000-iteration seed7/no-restart/no-
preconditioning repair; clean proposals are untouched. Held-out coordinates first
enter final scoring, never donor selection, context choice or refinement.

Report all per-seed raw/restored/final findings, heavy/backbone RMSD, movement,
source custody and lineage. This is one protein with multiple masks, not a family-
disjoint benchmark or qualified topology/trajectory guarantee. Validation seeds
remain unused until a separate frozen policy and independent cohort exist.

Implementation: `scripts/mac_conditioning_contracts.py`,
`scripts/preflight_mac_conditioning.py`, `scripts/run_mac_conditioning_study.py`.
Evidence: `tmp/mac-reconstruction-comparison-2026-11-11/conditioning-acd-attempt01/`.
