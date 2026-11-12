# THR2: local H mechanism and bounded matched-start transplants

## Scope and saved-only localization

Continues the [fresh native capture](thr2-native-capture.md) after the user's
instruction to continue research autonomously. Production code, force fields,
minimization settings and safeguards remain unchanged. Documentation scope is
this protocol/results note, a link from the capture note and the canonical
[handoff](../plans/pdbfixer-reconstruction-handoff.md); validate evidence and links
before publication. Original capture enrollment is finished, not extended.

`audit_thr2_angle_terms.py` analytically evaluates every native harmonic angle at
**all three retained whole-phase-1 starts**, with independent Reference force
readback. It does not optimize, place H or mutate forces (measurement groups only).
Finite-difference tests validate the analytic Cartesian derivative. Multi-H
methyl angles remain joint terms, never attributed twice.

| Start | Heavy-only angle Cβ projection | Target HB angle contribution | All angle terms |
|---|---:|---:|---:|
| unseeded | +909.351 | **−2,728.961** | −1,414.324 |
| seed 1 | +909.351 | **+2,978.857** | +4,469.238 |
| seed 2 | +909.351 | **+2,980.424** | +4,456.165 |

Units: kJ/mol/nm projected along the increasing-Cβ-volume Cartesian unit gradient.
Maximum analytic/native Cartesian-force discrepancy is **4.1×10⁻¹² kJ/mol/nm**.
The [complete angle audit](thr2-angle-terms-evidence.json) includes every focal
angle and H grouping. HB is a strong initial-force candidate, not yet proof of
whole-trajectory causality. Other H terms also differ.

## Prospective H-swap protocol

The [frozen plan](thr2-hswap-plan.json) is a byte-identical public metadata copy
of `tmp/pdbfixer-research/thr2-hswap-plan-01.json`; native payloads remain private.

The HB mask is chosen **after saved-only localization, before any transplant
minimization**. This is mechanistic follow-up on one poor-heavy THR2 context, not
an independently selected validation cohort. No new seed or native H placement.

Sources: the three retained original whole starts (unseeded, seed 1, seed 2),
with complete matching atom/bond identity, byte-identical Systems and identical
327 heavy coordinates. Transplant coordinates by exact atom identity only:

- **HB:** just `A/2/(empty insertion code)/THR/HB` (one H).
- **other-H:** every H except that HB (314 H); complementary partition.
- **unchanged:** each of the three starts without modification.

Enroll **all six ordered donor→recipient pairs × both masks**, plus three
unchanged baselines: **15 paths**, each once on captured-checkpoint OpenCL and
Reference XML routing, **30 paths / 60 native full-phase calls** total. No
outcome filtering, additional masks, donors, seeds, repeats or platform search.
Full-H substitution closes exactly to the donor start and needs no extra run.

Both native phases use **10 kJ/mol/nm / 1,000 iterations**, no solvent. Phase 2
chains from that path's own phase-1 return, not the original native phase-2
coordinates. Apply only the exact captured native restraint divide-ten change,
keeping anchors fixed; require serialized phase-2 System equality. No force
addition, H repair, reflection, coordinate restoration or production output.

Every actual starting coordinate must equal its frozen transplant exactly; do
not widen a transport tolerance. Retain starts/callbacks/returns privately.
Reporter copies native coordinates/gradients/stats and never requests early exit.
All outcomes, including D geometry, clashes, failures and nonexecutions, remain
reported rather than repaired or silently discarded.

Bounds: **4,096 callbacks per phase**, **1,800 s per path**, **14,400 s total**.
No retries or extensions. Independently audit exact masks/phase continuity,
System equality, stereochemistry and physical gradient/objective readback.
Unchanged Reference phase-1 controls must match the earlier saved-start replay
exactly; OpenCL endpoint variability is measured, not hidden.

Interpretation: successful bidirectional HB flips would support HB dependence in
these retained matched starts, not a universal mechanism, seed success rate,
production initializer or recovered November 6 trajectory. Complementary masks
test whether that interpretation survives the rest of H geometry. No outcome is
a production-qualification gate.

## Results: HB determines the branch within this finite matched-start set

Protocol checkpoint **`de12422`** precedes the transplant executions. All **30
paths / 60 native phases** completed and were independently audited. Initial
mask ownership is exact in every case; the 60 serialized phase Systems match
the original native phase parameters. Across **4,097 callbacks**, independent
same-platform mobile-gradient/objective errors are zero. No negative Cα or
returned/last-reported mismatch occurred. All 60 returns satisfy the native
gradient criterion; that does not qualify the structures as correct.

For both OpenCL and Reference, the final Cβ CIP matrix is:

| HB coordinates from ↓ / other 314 H from → | unseeded | seed 1 | seed 2 |
|---|---|---|---|
| unseeded | **S** | **S** | **S** |
| seed 1 | **R** | **R** | **R** |
| seed 2 | **R** | **R** | **R** |

Every phase-1 classification persists through phase 2. Swapping the single HB
rescues the unseeded S outcome with either R donor, and transfers the S outcome
to either R recipient. Transplanting the complementary 314 H preserves the
recipient HB's outcome. Initial total Cβ force projection changes accordingly:
unseeded with seed-1/2 HB is **+3,642.870 / +3,645.338 kJ/mol/nm** (Reference),
while seed-1/2 with unseeded HB is **−524.383 / −600.576**.

**Do not count reciprocal masks as independent samples.** HB(A←B) and
other-H(B←A) produce the exact same hybrid. The 30 numerical paths contain only
**nine distinct initial coordinate sets**, the complete 3×3 HB/remainder grid.
Six reciprocal Reference duplicate groups have exact phase-1 and phase-2
endpoints; all three unchanged Reference phase-1 controls exactly match the
earlier replay. OpenCL duplicate endpoints differ by up to **0.095902 Å**, without
changing any stereo outcome. The experiment provides finite-context causal
evidence for HB-coordinate dependence, not a success-rate estimate or broad
qualification. HB transplantation changes both its direction and bond length;
it does not isolate one harmonic-angle term from every other affected force.

[Complete transplant audit](thr2-hswap-evidence.json) retains every declared path,
input-equivalence group, sign crossing, terminal comparison and control. Private
coordinate/gradient traces remain under `tmp/pdbfixer-research/thr2-hswap-run-01/`.
Expanded independent private/public audits reproduce byte-identically; the first
audit remains retained rather than overwritten.

## Earlier localization: the saved initializer already chooses the conflicting side

The saved-only [H-initialization audit](thr2-h-initialization-evidence.json)
reconstructs **all 315 initial H coordinates / 945 saved Python random draws per
capture**. All three final Python RNG states match exactly. Maximum formula/native
position differences are **1.19–1.28×10⁻⁶ Å**, the measured OpenCL coordinate
transport scale, not double/single-precision bit identity. No native H placement
or minimization was rerun for this reconstruction.

The bound installed `Modeller.addHydrogens` formula initializes a bonded H from
the normalized sum of parent-minus-heavy-neighbor vectors, adds a positive-octant
`0.05*random()` perturbation, then normalizes again to 0.1 nm. For this retained
poor-heavy THR2 frame, the unperturbed HB formula gives four-ligand volume
**+2.547135 Å³**. The actual saved unseeded perturbation changes it to
**−1.571369 Å³ before H minimization**. Seed 1/2 remain positive
**+2.056973 / +0.686312 Å³**.

Two distinct geometric observables are reported, not conflated:

- Heavy-centre volume: `det(CA−CB, OG1−CB, CG2−CB)`.
- Four-ligand volume: `det(CA−HB, OG1−HB, CG2−HB)`.

In the unseeded capture, the heavy-centre volume remains **+0.198389 Å³** during
H-only minimization, while the four-ligand volume becomes **−5.789558 Å³**. That
HB geometry passes unchanged through Modeller return and the existing H-repair
bridge into whole minimization. Whole phase 1 then changes the heavy-centre
volume to **−2.472333 Å³**, while the four-ligand volume stays negative
**−8.148594 Å³**. The seeded controls have no such sign conflict.

This links the fresh event to a captured initializer-side choice, its persistence
through H-only minimization, and HB-dependent heavy inversion in matched-start
interventions. It is not evidence that removing all random perturbations is a
qualified initializer, that heavy reconstruction was sound, or that the lost
original November 6 run followed this exact trajectory. Runtime-declared source
provenance is not binary attestation. Existing THR24 terminal-leakage evidence
remains a separate mechanism.

## Validation and the current access blocker

The full fast run is **not green**: **3,082 passed / one failed**, three optional
Torch/Biotite skips and 33 slow deselections. The unchanged
`test_prepare_repairs_missing_hg_when_oxt_present` observed **OG–HG = 1.245813 Å**
against its existing **<1.2 Å** assertion. No threshold, production code or skip
was changed, and the test was not repeated until it passed.
The final focused research suite passed **26 tests**, including two saved-RNG
helper tests added after full-suite collection. Ruff, generated CLI/GUI specs,
agent-documentation checks, local links and public-payload checks passed.
Angle, transplant and initializer reports independently reproduce byte-identically.

The current source's distance-only O–H guard is `0.97*1.6 = 1.552 Å`, so that
guard alone would not flag this length. Other findings and the actual failed
structure have **not** been inspected: the read tool rejected access to the
external pytest directory. Do not invent a complete cause from the log or infer
that the HB study caused a regression in unchanged preparation code.

[Validation record](thr2-hswap-validation.json) binds the retained project-local
log/XML and unchanged source hashes. The original failing `hg_missing.pdb`,
`repaired.pdb` and `.dat` remain outside the allowed project scope at:

```text
/private/var/folders/rj/5v2182513dx99_nlm140sfl00000gn/T/pytest-of-lilpoundcake/pytest-1140/test_prepare_repairs_missing_h0/work/
```

**User input is now needed:** grant read/copy access to that specific directory
so the original failed artifacts can be preserved inside the project and audited
without reproducing or replacing the failure. No sandbox bypass was attempted.

**Subsequent continuation:** the user granted permission, but tool access remained
blocked and the user's copy command reported that the source path was absent.
Original raw artifacts remain unavailable. The separately bounded
[SER-HG fresh-capture stage](ser-hg-capture.md) preserves that loss and the original
failed validation result; it is not a recovered original run or retry-until-green.
