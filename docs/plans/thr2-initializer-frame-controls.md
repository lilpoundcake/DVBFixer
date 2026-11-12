# Prospective saved-RNG initializer frame controls

## Scope and documentation plan

After the user's autonomous-continuation instruction, return to the unresolved
THR2 mechanism on `research/thr2-initializer-frame-controls`, from `c6d9b21`.
O–H and GLYCAM fallback implementations/evidence stay unchanged and unmerged.
This is research-only source/formula analysis, not a new initializer capability,
FF change, native trajectory, public reconstruction route or upstream issue.

Documentation scope: this prospective contract, frozen JSON plan, scalar result
and validation records, a short research note and canonical handoff update.
Preserve all old hashes/outcomes and validate links/public payloads; do not link
to future files before creation. Private arrays remain in the already ignored
`tmp/pdbfixer-research/` tree. Pending future links: none.

## Source-supported question

The [three retained captures and HB transplants](../research/thr2-hswap.md)
already locate a four-ligand sign conflict in the unseeded HB initializer, before
H minimization. All 315 H/945 saved draws per capture were reconstructed. HB
transplantation changes both direction and radius; it does not isolate a single
angle term or qualify deterministic initialization.

The bound formula adds a vector with three positive `0.05*random()` components
in laboratory axes to a covariant heavy-neighbor-derived 0.1 nm direction.
Ask whether rotating an unchanged molecule while replaying exactly the same
saved draw triples can alter the HB side choice, and whether covariantly rotating
the perturbation eliminates that mathematical frame dependence. No endpoint or
chirality-rescue prediction is treated as a native trajectory result.

## Frozen finite enrollment (before evaluating transformed formulas)

- All three original labels: unseeded, seed 1, seed 2. Same original source
  coordinates, full atom/bond identity and saved RNG states as the old audit.
- **All 24 signed-permutation matrices with determinant +1**, including identity;
  no improper reflection or outcome-selected orientation.
- Two laboratory translations in nm: `(0,0,0)` and `(1.25,-0.75,0.5)`.
- **144 frames**, each evaluated with all three formula arms for all 315 H:
  **432 formula-arm records / 136,080 initialized-H evaluations**.
  Additionally perform **945 identity-formula checks** through the already
  audited helper; these are bookkeeping controls, not additional frame arms.
  These are correlated deterministic transformations of one structure/three
  saved draw sets, not 432 physical samples, native runs or independent trials.
- Arms: (1) native fixed-lab saved perturbation, (2) same saved vector rotated
  with the molecule (explicit counterfactual covariance control), (3) zero
  perturbation (mathematical control, not a qualified replacement initializer).

Replay each saved draw stream exactly, with exact final RNG-state equality and
no new seeds/global RNG mutation. Use the already audited heavy input coordinates
and bonded-neighbor order. Native-formula identity-frame positions must equal
the existing audited formula exactly; compare with saved OpenCL coordinates as
the already measured transport discrepancy, not double/single bit identity.

## Measurements and engineering assertions

Primary observable remains the distinct volumes
`det(CA-CB, OG1-CB, CG2-CB)` and `det(CA-HB, OG1-HB, CG2-HB)`.
Their sign conflict is **not a CIP assignment** or proof of a heavy inversion.
Heavy geometry is only rigidly transformed and must not change chirality.

Secondary controls enroll, by graph before outcomes, every carbon with exactly
three heavy neighbors and exactly one attached H. Retain all their four-ligand
sign conflicts, not just THR2. Multi-H methyl groups are initialized and included
in covariance errors but are not single-H tetrahedral sign tests.

Inverse-transform every proposal. Report full-H and target-H discrepancies from
the identity formula. Covariantly rotated-noise and zero-noise controls must
agree within **1e-12 nm**; translation-only paired comparisons use the same bound.
This tolerance concerns deterministic double-precision algebra, not native GPU
transport or scientific acceptance. Degenerate raw neighbor sums ≤1e-12 nm
refuse explicitly. Volumes within **1e-10 Å³** of zero are reported degenerate,
not classified positive/negative. Independently audit determinants and formula
results (determinant agreement ≤1e-10 Å³); retain any assertion failure without
changing tolerances/orientations.

Keep every frame/arm and refusal. Counts over the 24 orientations are not a
probability estimate for uniformly random molecular orientations or fresh RNG
sampling. Do not search for a favorable seed, angle, perturbation strength or
new generator. No Context, energy/force query, optimizer, native H placement,
repair, reflection, table lookup or candidate PDB publication is permitted.

Pathwise fixed-draw covariance and invariance of a sampling distribution are
different questions: even isotropic laboratory noise need not give identical
rotated trajectories for one fixed draw vector. The source's positive cube is
not an isotropic distribution, but these three saved streams cannot estimate
its orientation-dependent inversion probability. Do not infer a general
initializer defect or a production remedy from pathwise noncovariance alone.

Bound execution to **120 s**, exactly the declared 144 frames/three arms and saved
draw counts. No retries after a failed scientific assertion, no added frames or
adaptive extension. Source/fixture/runtime provenance and original artifacts
remain hash-bound. Future initializer intervention or minimization requires its
own separately frozen protocol; this stage does not authorize production use.
