# Apple Silicon Protpardelle 100-Case Soak

## Decision

The frozen first-100 operational soak passed: 100/100 cases completed native
MPS inference, generated-only OpenMM CPU refinement, validation, and atomic
publication without fallback, process failure, memory-pressure termination, or
temporary-directory leakage. Scientific failures remain in the denominator;
validation-first selection passed 91/100 cases.

The user subsequently authorized extension to the full frozen 231-case cohort.
That extension is independently frozen in
[`small-diffusion-apple-full-231.json`](small-diffusion-apple-full-231.json).

## Frozen Inputs

- Manifest: [`small-diffusion-apple-100-soak.json`](small-diffusion-apple-100-soak.json).
- Membership: the first 100 cases in `small-diffusion-full-followup-v2.json`
  manifest order, committed before inference in `355717d`.
- Protpardelle revision: `ee378400f25b801fa481028000f9060183d7fb4c`.
- PyTorch 2.6.0 float32 MPS, 500 steps, fallback disabled, one fresh process
  per case.
- Refinement: OpenMM CPU, generated-only, eight restarts, no resampling.
- Host inventory SHA-256:
  `a2ea71ac6fa5833cfe15bd98e9519188c911834dd934b2d72331af57961721e4`.
- Machine-readable ignored result:
  `.artifacts/apple-diffusion/apple-100-result.json`, SHA-256
  `45bd3b38cb88aaad1d2585a15633155d0bed7c8f129d57d1b423da54e029bbe8`.

## Results

- Operational completion: **100/100**.
- Raw MPS hard-gate passes: 42/100.
- Refined passes: 91/100.
- Validation-first selected passes: **91/100** (48 gap-5, 43 gap-10).
- Median selected gap-backbone RMSD: 0.534 A.
- Failed selected cases: `9g9t`, `9e51`, `9ges`, `10rx`, `28op`, `22uy`,
  `24rc`, `13io`, and `21vq`.

Every published candidate retained exact fixed-heavy coordinates. Scientific
hard-gate failures did not stop or replace later cases.

## Runtime And Memory

- Median MPS inference wall time: 18.21 s.
- Median CPU refinement wall time: 12.33 s.
- Linear-interpolation p95 raw-plus-refinement time: 103.14 s.
- Maximum raw-plus-refinement time: 464.37 s.
- Peak process RSS: 912,080,896 bytes (4.72% of physical unified memory).
- Peak sampled MPS tensor allocation: 289,383,168 bytes.
- Peak sampled Metal driver allocation: 355,254,272 bytes.
- Residual staging directories after completion: zero.

The memory counters are overlapping views of unified memory and are not summed.
The long tail comes from CPU refinement rather than MPS inference.

## Structure Export

Input, experimental reference, and validation-first output structures are
available under:

```text
.artifacts/apple-diffusion/full-231-structures/<case-id>/
  <case-id>_input.pdb
  <case-id>_reference.pdb
  <case-id>_output.pdb
  <case-id>_metadata.json
```

Failed cases retain the refined failure artifact as `<case-id>_output.pdb` and
declare `selection: failed-refined` plus the hard-gate failures in metadata.
