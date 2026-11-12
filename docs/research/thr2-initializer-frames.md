# THR2 initializer: fixed-draw coordinate-frame controls

## Outcome

Prospective checkpoint **`2b58d4e`** precedes all transformed source evaluations.
The [frozen protocol](../plans/thr2-initializer-frame-controls.md) and
[JSON plan](thr2-initializer-frames-plan.json) retain all three saved RNG streams,
24 proper signed-permutation rotations, two translations and three formula arms.
All **144 frames / 432 formula-arm records** completed, covering all **315 H per
frame/arm**. No Context, optimizer, native H-placement call, new random stream,
force-field intervention, geometry repair or PDB publication was used.

[Independent complete audit](thr2-initializer-frames-evidence.json) reports all
outcomes and arithmetic controls. These are correlated transformations of one
poor-heavy THR2 context, **not 432 physical samples or native trajectories**.
The [validation record](thr2-initializer-frames-validation.json) binds retained
source/evidence digests and full-suite results.

At zero translation, the target heavy-centre/four-ligand sign conflicts are:

| Saved draw stream | Fixed laboratory perturbation | Covariantly rotated perturbation | Zero perturbation |
|---|---:|---:|---:|
| unseeded | **6/24** | **24/24** | 0/24 |
| seed 1 | 0/24 | 0/24 | 0/24 |
| seed 2 | **3/24** | 0/24 | 0/24 |

Translation does not change these classifications. The 56 graph-selected
single-H carbon centers per source were all retained as secondary controls;
**only THR2/HB** has a sign conflict in any arm/frame. All these centers remain
nondegenerate under the predeclared 1e-10 Å³ rule. Proper rotations preserve
heavy geometry: no heavy chirality inversion occurs in this algebra experiment.

## What the result means

In these saved draws, changing only the coordinate frame can alter the native
initializer's HB side choice even though the molecule's heavy geometry is the
same. The bound source adds a positive-octant perturbation in laboratory axes
to the heavy-neighbor-derived direction. A proper rigid frame change does not
rotate that fixed saved laboratory vector.

The covariant-offset control does **not rescue** the unseeded initializer: it
correctly preserves its original sign conflict under every rotation. It removes
pathwise frame dependence while retaining that original bad choice. Zero-offset
HB is sign-compatible throughout these frames, but was **not minimized** and is
not a qualified initializer or a demonstrated remedy for native Cβ inversion.

Keep the observables distinct:

- Heavy-centre volume: `det(CA−CB, OG1−CB, CG2−CB)`.
- Four-ligand volume: `det(CA−HB, OG1−HB, CG2−HB)`.
- Their opposite signs indicate this geometric conflict, **not a CIP label**.
  This stage does not measure a dynamic R→S transition.

Also distinguish **pathwise covariance at a fixed saved draw** from **invariance
of a sampling distribution**. Even isotropic laboratory noise need not yield a
pathwise-rotated result for the same fixed draw vector. The source's positive
cube is nonisotropic, but three saved streams and this finite rotation group do
not estimate its orientation-dependent inversion probability. Counts 6/24 and
3/24 are neither failure rates nor evidence that deterministic H placement is
generally more accurate. No favorable frame was selected for minimization.

### Why the global zero-offset arm is not a physical initializer

A separately labelled [retrospective readback](thr2-initializer-zero-offset-limit.json)
of the already retained arrays finds **139 exactly coincident same-parent H
pairs in every zero-offset frame**. Sibling H share the same heavy-neighbor
direction when their perturbations are all removed. No new proposal, frame,
formula evaluation or native call was added to obtain that diagnostic; it is
not a predeclared acceptance criterion retrofitted to the experiment.

The noisy arms have no exactly coincident sibling pairs in these arrays, but
some are very close (minimum fixed-lab pair distance 0.000123 Å), so even those
initial arrays are not certified physical geometries. The algebra pass means
only the declared coordinate/formula controls agree. Removing noise globally
is **not a viable production conclusion** from this study; a future HB-only
counterfactual must retain the other 314 H and their saved perturbations.

## Independent validation

The runner uses the saved 945-draw stream per original capture and exactly
matches its final Python RNG state. The 945 identity-helper checks use the
previously audited source formula; its measured formula/OpenCL discrepancy stays
at the historical ~1.2e-6 Å transport scale, not double/single bit identity.

The independent NumPy audit re-evaluates all transformed formulas/determinants,
all proper rotations and all translation pairs. Maximum observed errors:

- Formula-coordinate agreement: **6.305e-16 nm**.
- Translation-pair inverse-coordinate agreement: **5.556e-15 nm**.

Both are below the predeclared **1e-12 nm** algebra bound. Covariant-offset and
zero-offset proposals agree after inverse transformation within that bound.
Native identity-frame formula checks are exact. No tolerance, frame, seed,
perturbation strength or degeneracy rule was changed after source evaluation.

Five toy engineering tests check the complete 24-rotation group, exact native
offset algebra, covariant/translated round trips, explicit fixed-lab
noncovariance and singular-neighbor refusal. They are not additional captured
sources or native samples. All private proposal arrays and complete execution
evidence remain under `tmp/pdbfixer-research/thr2-initializer-frames-run-01/`.
Independent offline audit readback is byte-identical. The complete fast suite
passed **3,207 tests**, with three optional Torch/Biotite skips and 33 slow
deselections; the five focused toy controls, Ruff, strict-target mypy, generated
CLI/GUI specs and agent-map checks also passed. Production and original SER
regression assertions are unchanged from `c6d9b21`.

## Remaining mechanism gate

The next separate question is dynamic: if **only HB's perturbation** is removed
before the original H-only minimization, does a sign-compatible initializer
remain compatible and change the ensuing heavy Cβ branch? That needs a new
frozen counterfactual protocol, unchanged native FF/minimizer settings, all
three saved sources, original controls, exact bridge/identity checks and retained
negative outcomes. This algebra result does not answer it.

A [saved-only admission inspection](thr2-hbridge-admission.json) confirms that the original historical H-repair
bridge adds **zero repairs** on each saved H return and exactly reproduces all
three original whole-phase-1 coordinate arrays. No new Context was created for
that inspection. Future counterfactual H returns still need the original
historical bridge and independent auditing; do not silently substitute the
newly changed production O–H guard into old captured physics.

Production remains unchanged from `c6d9b21`, with both approved fixes still on
separate unmerged branches. The original THR24 terminal leakage, lost November 6
THR2 trajectory, lost original SER endpoint and empirical-table production and
license gates remain separate and unchanged.
