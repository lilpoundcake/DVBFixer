# Apple Silicon Protpardelle Full 231-Case Results

## Decision

The user-authorized full frozen cohort completed on Apple Silicon: **231/231
cases reached a terminal, provenance-recorded outcome**. Native MPS inference
completed for every case with PyTorch CPU fallback disabled. Validation-first
selection passed 213/231 cases (92.2%): 121/128 gap-5 and 92/103 gap-10.

This establishes operational feasibility on an 18 GB Apple M3 Pro. It does not
show that MPS is faster than CPU: the short-case CPU reference remained about
three times faster, and the runtime long tail came from OpenMM CPU refinement.
The exact checkpoint size, model parameter count, architecture summary, and
cross-engine comparison are recorded in
[`diffusion-model-size-and-run-parameters.md`](diffusion-model-size-and-run-parameters.md).

## Frozen Protocol

- Manifest: [`small-diffusion-apple-full-231.json`](small-diffusion-apple-full-231.json).
- Membership: all 231 cases from `small-diffusion-full-followup-v2.json`, in
  frozen manifest order; 128 gap-5 and 103 gap-10 cases, with no replacement.
- Protpardelle revision: `ee378400f25b801fa481028000f9060183d7fb4c`.
- PyTorch 2.6.0 float32 MPS, 500 steps, fallback disabled, one fresh inference
  process per case.
- Generated-only OpenMM CPU refinement, eight restarts, no diffusion resampling.
- Validation-first selection: valid refined candidate, otherwise valid raw
  candidate, otherwise a retained failure.
- Host inventory SHA-256:
  `a2ea71ac6fa5833cfe15bd98e9519188c911834dd934b2d72331af57961721e4`.

The ignored machine-readable aggregate is
`.artifacts/apple-diffusion/apple-full-231-result.json`, SHA-256
`5e54d25e7830eb5842a20aa5cb647244ff0d76ff66ac0553caddd0a5ab30200c`.

## Scientific Results

| Stage | All | Gap 5 | Gap 10 | Median selected backbone RMSD | Median selected all-heavy RMSD |
|---|---:|---:|---:|---:|---:|
| Raw MPS | 96/231 | 67/128 | 29/103 | — | — |
| Refined CPU | 212/231 | 120/128 | 92/103 | — | — |
| Validation-first selected | **213/231** | **121/128** | **92/103** | **0.506 Å** | **1.356 Å** |

One case (`35zv`) retained its valid raw candidate after refinement failed
validation. Eighteen selected cases failed: `9g9t`, `9e51`, `9ges`, `10rx`,
`28op`, `22uy`, `24rc`, `13io`, `21vq`, `13de`, `9eam`, `9ge5`, `12zn`,
`9ejk`, `11ex`, `9gyi`, `9ddr`, and `9jjb`.

The corresponding Linux CUDA follow-up reported 101/231 raw and 217/231
selected passes. Apple MPS therefore produced five fewer raw passes and four
fewer selected passes on the same frozen membership. Exact training membership
remains unresolved, so both results remain descriptive rather than backend
selection evidence.

### `9ge5` Chirality Failure

The first `9ge5` CPU refinement raised `ChiralityError` for `H/ARG223` and
published no candidate, preserving the zero-D-Cα invariant. A diagnostic rerun
of the same nominal refinement completed, demonstrating the already documented
OpenMM CPU nondeterminism, but remained invalid and was excluded from selection.
The official outcome is the first scientific refinement failure; the raw MPS
candidate also failed hard gates. The export therefore labels this case
`failed-raw` and records the refinement failure in metadata.

## Operational Results

- Terminal outcomes: **231/231**.
- Median MPS inference wall time: 17.87 s.
- Median CPU refinement wall time: 9.61 s (230 timed published refinements).
- Linear-interpolation p95 raw-plus-refinement time: 63.02 s.
- Maximum raw-plus-refinement time: 464.37 s.
- Peak process RSS: 912,080,896 bytes (4.72% of physical unified memory).
- Peak sampled MPS tensor allocation: 289,383,168 bytes.
- Peak sampled Metal driver allocation: 355,254,272 bytes.

The memory counters are overlapping unified-memory views and are not summed.
No MPS fallback, wrong-device tensor, memory-pressure termination, identity
mapping failure, incomplete atom set, or non-exact published fixed atom was
observed.

## Structure Export

All 231 cases are exported under:

```text
.artifacts/apple-diffusion/full-231-structures/<case-id>/
  <case-id>_input.pdb
  <case-id>_reference.pdb
  <case-id>_output.pdb
  <case-id>_metadata.json
```

There are exactly 231 input, 231 reference, and 231 output PDB files. Failed
cases remain clearly marked in metadata and are not presented as valid selected
structures. The export inventory is `export-manifest.json`, SHA-256
`1a246f34469867332b948ffe2456a3ecc0f3215e7222552ba7358f3beaacd738`.
