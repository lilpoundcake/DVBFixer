# THR2: HB-only initialization changes the native Cβ branch

## Outcome and scope

The [prospective protocol](../plans/thr2-hb-only-initializer.md) and
[frozen plan](thr2-hb-only-plan.json) were committed at **`e822dbc` before any
intervention execution**. All **12 paths / 36 native minimizer calls** completed
once, with no new seeds, replacements, retries or outcome-selected extensions.
The [complete independent audit](thr2-hb-only-evidence.json) retains every path,
callback measurement, control and phase transition.
The [validation record](thr2-hb-only-validation.json) binds source/evidence hashes,
native stopping limitations and the final full-suite result.

Only the initial coordinate of `A/THR2/HB` is changed in the intervention arm.
The other **314 initialized H and all 327 heavy coordinates** are initially
copied exactly; all original H remain mobile in the H-only minimization, while
the heavy atoms remain pinned. Both whole-system phases then use their original
mobility and exact original System parameters. This is a saved-state coordinate
counterfactual, not a production `Modeller.addHydrogens` implementation.

Final Cβ CIP outcomes in both captured OpenCL and Reference modes:

| Saved start | Unchanged initialization | Only HB made unperturbed |
|---|---|---|
| unseeded | **S** | **R** |
| seed 1 | R | R |
| seed 2 | R | R |

The unchanged unseeded path reproduces the fresh captured R→S Cβ inversion in
whole phase 1, persisting through phase 2. Changing only its initial HB prevents
that inversion under these exact settings. Seeded controls retain R in either
arm. The two numerical modes are not independent biological samples: this is
**one poor-heavy context / three saved H starts / six proposals**, not a success
rate estimate over twelve structures or newly sampled seeds.

## The causal chain that is actually supported

The preceding [frame controls](thr2-initializer-frames.md) establish only formula
behavior, and the earlier [HB transplants](thr2-hswap.md) localize dependence on
the post-H HB coordinates within a finite donor grid. This new experiment moves
the intervention upstream, before the original H-only minimizer, while retaining
all other initialized coordinates and the original native physics.

For the unseeded start, the Reference path illustrates the sequence:

| Stage | Unchanged HB four-ligand volume | HB-only-unperturbed four-ligand volume | Heavy Cβ volume |
|---|---:|---:|---|
| Before H-only | −1.571369 | +2.547135 | +0.198389 in both |
| After H-only | −5.787773 | +6.891830 | +0.198389 in both |
| After whole phase 1 | −8.143315 | +8.036684 | −2.469279 vs +2.454817 |
| After whole phase 2 | −8.137384 | +8.067272 | −2.446946 vs +2.457067 |

Volumes are Å³. Heavy-centre and four-ligand volumes are different observables;
their signs are not substituted for RDKit CIP assignments. Cα remains S at the
target and no negative Cα volume is observed anywhere in any recorded path.

The corresponding OpenCL intervention gives the same R outcome, with numerical
trajectory differences retained. The newly computed OpenCL HB point is slightly
quantized (initial four-ligand volume +2.547136 rather than +2.547135); this lies
inside the prospectively declared requested-coordinate rounding bound. Do not
assert bit-identical OpenCL endpoints or replace the Reference values with them.

This supports **initializer-coordinate → H relaxation branch → heavy Cβ branch**
dependence in the retained context. It does not isolate one harmonic angle from
all other affected force terms. The initial point is changed, other H can respond
during native relaxation, and the full FF acts afterward. It does not prove a
universal FF defect, demonstrate a generally correct H generator, establish
rotamer accuracy or reconstruct the lost November 6 trajectory.

## Controls, physics and independent readback

- All **36 phase System serializations match the original native Systems**.
  H-only uses the original **1.0 kJ/mol/nm / 50 iterations**; whole phases use
  **10 kJ/mol/nm / 1,000 iterations** each, with no solvent. Phase 2 uses only
  the exact native same-anchor divide-ten restraint transition.
- Exact requested HB masks and unchanged initial copied coordinates are verified
  for all paths. All heavy coordinates remain exact throughout H-only callbacks
  and return. The largest measured requested/applied component discrepancy is
  **4.25814e-8 nm**, within the predeclared OpenCL float32 rounding bound.
- All **three unchanged Reference H-only endpoints exactly equal the previous
  retained Reference saved-call controls**. Whole phases chain from their own
  current H returns; no old native return is spliced into a counterfactual.
- The historical H-repair bridge is independently replayed on every H return:
  **zero repairs in all twelve paths**, exact coordinates and heavy identity.
  The newly implemented production O–H guard is not substituted into old physics.
- Across **2,174 callbacks**, independent same-platform physical mobile-gradient
  and objective readbacks have **zero maximum error**. No returned/last-callback
  position mismatch occurs in any phase. No source forces, constraints, walls,
  optimizer settings or anchors are tuned; no emergency reflection is applied.
- The complete independent offline audit repeats byte-identically. This is a
  readback of already executed paths, not an extra minimization or native sample.

### Stopping and geometry limitations stay visible

All twelve H-only calls produce 50 callbacks and **none satisfies the measured
native gradient criterion at return**. We did not extend that original budget
or tune a stopping rule after observing the outcome. All 24 whole-system returns
satisfy the native gradient criterion; this is not a global-minimum certificate.

After H-only relaxation, maximum O–H distances across these returns remain
below 0.982 Å and minimum same-parent H–H distances exceed 1.636 Å. These are
descriptive metrics, not comprehensive steric/protonation qualification. Raw
coordinates retain other contacts and intermediate poor geometry for inspection.

The prior global zero-offset mathematical control is **not this intervention**:
it gives 139 exact same-parent H coincidences per frame. Here only one single-H
parent's point changes, while the other 314 noisy initialized H are retained.
Do not turn this result into blanket removal of random perturbations.

## Engineering and publication boundary

Six new mask/transport tests (eleven with the preceding frame toy tests) passed.
They cover the exact single-H mask, immutable input arrays, unchanged controls,
undeclared-arm refusal, exact Reference transport, bounded OpenCL rounding and
nonfinite-coordinate refusal. They do not execute experimental minimization.
Original phases, aliases and bond graphs were admitted against hash-bound old
captures and the canonical shared variant map before native execution.

The complete fast suite passed **3,213 tests**, with three optional
Torch/Biotite skips and 33 slow deselections. Ruff, strict-target mypy, generated
CLI/GUI specs, agent-map and public-payload/link checks passed. These engineering
passes do not override the nonconverged H-only returns or establish a broadly
qualified replacement initializer. Runtime provenance is declared/hash-bound;
it is not an independent attestation of GPU binary identity.

All private source coordinates, requested/applied arrays, Systems, native states,
callback gradients, partial-path infrastructure and logs remain under
`tmp/pdbfixer-research/thr2-hb-only-run-01/`. Tracked records contain identities,
hashes and scalar measurements only; no candidate PDB is published.

Production remains unchanged from `c6d9b21`, and `main` is untouched. An integrated
H-initialization guard would be a **new behavior-changing production design**,
not automatic promotion of this research script. It requires separate approval,
source/mask admission and broader qualification on independent contexts,
especially observed-H preservation, well-formed heavy geometry, near-planar
centers, other chemistries, actual native FF routes and stopping behavior. The
underlying poor-heavy reconstruction and empirical-library production/license
gates remain separate unresolved work.
