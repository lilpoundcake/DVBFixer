# Draft technical appendix — not an issue-ready human summary

**Publication gate:** do not submit this file unchanged. OpenMM's
[AI policy](https://github.com/openmm/openmm/blob/3e34aa2091bfff5306206b4239c4654776342fd8/AI_POLICY.md)
requires an opening summary written by the human submitter in their own words,
and their review of any AI-generated material. That opening has not been
provided. Everything below is explicitly **AI-generated technical material**,
prepared from locally executed experiments; it must follow the human opening.
Do not claim a human review or answer maintainer questions automatically.

Suggested issue title (AI-generated, subject to human review):

**Constraint penalty gradient appears reversed in GPU minimizer CPU fallback
(three-particle OpenCL reproducer)**

---

## AI-generated technical appendix

### Observed behavior and expected derivative

A three-particle synthetic System produces a public MinimizationReporter
gradient inconsistent with the derivative of its objective (System energy plus
positive constraint penalty) on the tested OpenCL single-precision runtime.
The source-sign hypothesis
from CommonMinimizeKernel::evaluateCpu, after float32 transport, matches all nine
reported components exactly. CPU/Reference and eleven other factorial controls
pass the physically correct derivative checks. No native source was patched.

For the constrained particle 0 x component in the successful trigger:

| Quantity | kJ/mol/nm |
|---|---:|
| Analytic physical + positive penalty gradient | 199.4642138481195 |
| Independent termwise central finite difference | 199.46421384833772 |
| Native reporter gradient | 200.52517700195312 |
| Analytic physical - penalty gradient | 200.52517652511037 |

The discrepancy~1.060963 exceeds the predeclared strict comparison bound~0.005011
PLUS the additional single-precision constraint-arithmetic allowance~0.762941.
Matching the opposite sign is not accepted as a physically correct derivative.
Physical and penalty energy errors at this node are both0.

### Tested environment

- macOS 26.6.2, arm64, Apple M3 Pro GPU.
- Python 3.12.13; OpenMM self-reports `8.5.2.dev-36a30cb`, revision
  `36a30cbca54e727b216b606f3c011b67201eb8b4`.
- Existing micromamba-managed Python environment; no installation or package
  modification performed for these tests. Package channel/build provenance
  beyond the runtime's version strings has not been independently attested.
- Available platforms: Reference, CPU, OpenCL; plugin-load-failure list empty.
- Affected test: OpenCL, `Precision=single`, `DeviceIndex=0`,
  `OpenCLPlatformIndex=0`, `UseCpuPme=false`, `DisablePmeStream=false`.
- Portable script requires only stdlib+OpenMM, no NumPy, biological inputs,
  force-field templates, dvbfixer installation or private files. Independent
  factorial auditing used NumPy for analytic/finite-difference checks.
- `python -m openmm.testInstallation` was **not run** as part of this limited
  protocol. The platform controls below and plugin metadata are not that command's
  output and must not be presented as an installation-test pass.
- Neither current upstream master nor other releases/hardware were tested.
  Runtime metadata is not source/binary identity attestation.

### Complete portable reproducer

[Standalone script at immutable commit 32be60e](https://github.com/lilpoundcake/DVBFixer/blob/32be60e/scripts/openmm_cpu_sign_minimal.py)
([raw download](https://raw.githubusercontent.com/lilpoundcake/DVBFixer/32be60e/scripts/openmm_cpu_sign_minimal.py)).
Save as `openmm_cpu_sign_minimal.py` in an OpenMM environment, then run:

```bash
python openmm_cpu_sign_minimal.py --run-native --platform OpenCL
```

Without `--run-native`, the script describes only and does not create a Context.
CPU/Reference controls are available via `--platform CPU` / `--platform Reference`.
The exact compact OpenCL entry was run once after the factorial study as an
explicit packaging check: exit0, returned=true, one callback, same scalar witness
above, source-sign-float32 maximum error0. It was not counted as another
independent scientific case or qualification seed. No changes or retries were
needed. The full study contains24 calls plus this one separately recorded entry
verification, not25 independent scientific cases.

The System consists of three 12-Da particles, a harmonic bond (r_eq=0.12 nm,
k=10000 kJ/mol/nm²), a pair distance constraint of 0.1 nm and an independent quadratic
stress force on particle 2. All masses are positive; there are no virtual sites
or nonbonded interactions. A VerletIntegrator is used without integration,
constraint tolerance 1e-5; the minimizer runs with tolerance 10, maxIterations=1.
This is a gradient diagnostic, not a convergence or global-minimum claim. The
code contains all synthetic constants and initial positions; no input files are
needed. The reported gradient magnitude1,201,164,672 is above detector2^30 but
below signed fixed-point capacity2^31, avoiding the separate excessively-large-
force conversion problem.

### Source hypothesis, not a directly traced internal event

With `delta=(positions[y]-positions[x])/r`, the positive penalty
`0.5*K*(r-r0)^2` has `grad_y += K*(r-r0)*delta` and
`grad_x -= K*(r-r0)*delta`.

The inspected
[exact-revision CPU fallback source](https://github.com/openmm/openmm/blob/36a30cbca54e727b216b606f3c011b67201eb8b4/platforms/common/src/CommonMinimizeKernel.cpp#L617-L639)
instead subtracts from y and adds to x, after adding positive penalty energy.
It later [uploads hostGrad](https://github.com/openmm/openmm/blob/36a30cbca54e727b216b606f3c011b67201eb8b4/platforms/common/src/CommonMinimizeKernel.cpp#L658-L663)
to the device gradient buffer. This explains the observed sign/float32 pattern,
but no internal CPU-fallback event was directly instrumented. The observed
public-reporter/objective inconsistency is firmer than that internal-path
attribution. No proposed patch has been tested or submitted.

### Controls and preserved negative experiments

The successful representable-force series had 12 cases, each one callback:
OpenCL/CPU/Reference × constraint on/off × below/above detector. The constrained
OpenCL trigger rejects the physical gradient while the source-sign prediction
passes strict component comparison; all 11 other controls pass physical
gradients, physical/penalty energies and termwise finite differences. The
unconstrained OpenCL trigger and constrained below-threshold OpenCL control
both pass, distinguishing the sign concern from general out-of-range force
conversion.

An earlier 12-case oversized-force series did NOT produce the intended sign
witness. Its OpenCL unconstrained/constrained high-force gradients were already
incorrect, and one low-force constrained callback failed the old energy bound.
Those negative controls are preserved, not discarded or relabeled as passes.
Constants and starting points for the subsequent representable-range study
were frozen separately before execution; no comparison tolerances were widened.

Original comparison envelopes were fixed in advance: gradient1e-3 plus2e-5
relative per component; physical energy1e-3 plus2e-6 relative; penalty energy
1e-3 plus2e-5 relative. OpenCL constraint endpoints additionally report the
fixed arithmetic allowance `32*eps32*abs(K)*(r+r0)`. The sign mismatch remains
outside it. Every successful-series callback passes independent termwise
central finite differences (h1e-6nm), avoiding cancellation between small pair
energy and large independent stress energy.

### Related upstream reports and limitations

- #5394 concerns unconstrained CUDA fallback performance; this report is a
  distinct correctness observation involving constraint gradients.
- #5434 / merged PR#5435 concern Apple M3 out-of-range float-to-long conversion.
  The successful trigger here stays inside signed fixed-point range. The older
  oversized negative series should not be filed as a duplicate of that known
  overflow report.
- No exact constraint-sign duplicate was found in the issue/PR searches performed;
  this is not an exhaustive guarantee that none exists.
- No claim of effects on realistic structures, causality for any other geometry
  defect, installed-binary/source identity, or every platform/release.
- No private structures, measured coordinate/RNG arrays, credentials, local paths
  or proprietary parameter tables are needed or included.

Please check whether the penalty-gradient sign in `evaluateCpu` should be
reversed to agree with the positive objective, and whether a constrained,
representable-force overflow regression would be appropriate. This is a bug
report request, not a claim that a dependency fix has already been validated.
