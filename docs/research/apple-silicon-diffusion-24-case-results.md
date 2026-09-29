# Apple Silicon Protpardelle 24-Case Results

## Decision

The preregistered Apple 24-case portability gate passes at its exact minimum:
21/24 validation-first selected candidates passed every hard gate. All 24 cases
completed MPS inference and CPU refinement operationally. Peak process RSS was
4.70% of physical unified memory and the linearly interpolated p95 complete
raw-plus-refinement time was 92.47 seconds, below the frozen 300-second bound.
No new hard-gate failure class appeared. The plan therefore permits the
100-case operational soak, while the result does not establish an MPS speed
advantage over native CPU.

## Frozen Protocol

- Manifest: [`small-diffusion-apple-24-pilot.json`](small-diffusion-apple-24-pilot.json).
- Membership: the previously frozen 12 gap-5 and 12 gap-10 cases, unchanged and
  committed before Apple inference in `ece8188`.
- Protpardelle revision: `ee378400f25b801fa481028000f9060183d7fb4c`.
- PyTorch 2.6.0 float32 on MPS; `PYTORCH_ENABLE_MPS_FALLBACK=0`.
- 500 steps, step scale 1.2, churn 200.0, request seed, no MPNN, no per-step
  reinjection.
- One fresh process per raw case. Generated-only OpenMM CPU refinement used
  eight restarts and did not resample the diffusion model.
- Selection: refined candidate when valid, otherwise raw candidate when valid;
  every other case remains a failure.

The ignored machine-readable result is
`.artifacts/apple-diffusion/apple-24-result.json`, SHA-256
`d584fdec9c0ee6ed4df4ff30ccb8eeda9f0a170aa3f179091be7cc2312da6cd1`.

## Results

| Case | Length / gap | Raw MPS | Refined CPU | Selected | MPS wall | Refine wall | Selected backbone RMSD |
|---|---:|---:|---:|---|---:|---:|---:|
| `9hcf` | 321 / 5 | pass | pass | refined | 19.14 s | 17.31 s | 0.435 A |
| `8x9z` | 94 / 5 | fail | pass | refined | 18.12 s | 3.06 s | 0.322 A |
| `9e85` | 159 / 5 | pass | pass | refined | 19.10 s | 3.38 s | 0.203 A |
| `30ie` | 415 / 5 | fail | pass | refined | 19.61 s | 43.24 s | 0.323 A |
| `11mi` | 176 / 5 | fail | pass | refined | 19.11 s | 5.20 s | 0.304 A |
| `24mc` | 461 / 5 | fail | pass | refined | 20.18 s | 110.55 s | 2.035 A |
| `9dvi` | 77 / 5 | fail | pass | refined | 18.89 s | 3.80 s | 1.813 A |
| `9jft` | 282 / 5 | fail | pass | refined | 19.33 s | 21.35 s | 0.640 A |
| `23vl` | 289 / 5 | fail | pass | refined | 19.51 s | 22.28 s | 0.582 A |
| `9ges` | 177 / 5 | fail | fail | fail | 19.07 s | 78.63 s | — |
| `45iy` | 389 / 5 | pass | pass | refined | 19.54 s | 18.11 s | 0.452 A |
| `28ye` | 99 / 5 | fail | pass | refined | 19.09 s | 3.04 s | 1.604 A |
| `9eko` | 452 / 10 | pass | pass | refined | 19.99 s | 34.02 s | 3.158 A |
| `9guc` | 82 / 10 | pass | pass | refined | 18.36 s | 2.82 s | 0.321 A |
| `9g9t` | 205 / 10 | fail | fail | fail | 19.23 s | 22.17 s | — |
| `9jq3` | 387 / 10 | fail | pass | refined | 19.64 s | 23.04 s | 0.477 A |
| `9hab` | 385 / 10 | fail | pass | refined | 19.73 s | 30.38 s | 0.495 A |
| `8v9v` | 167 / 10 | pass | pass | refined | 19.13 s | 5.99 s | 0.392 A |
| `8uui` | 295 / 10 | fail | pass | refined | 19.54 s | 17.10 s | 1.401 A |
| `11fs` | 112 / 10 | fail | pass | refined | 18.80 s | 3.43 s | 0.945 A |
| `9e51` | 215 / 10 | fail | fail | fail | 19.11 s | 14.96 s | — |
| `9g1v` | 19 / 10 | fail | pass | refined | 19.84 s | 0.69 s | 2.628 A |
| `9e2t` | 201 / 10 | pass | pass | refined | 19.11 s | 8.34 s | 1.451 A |
| `9cf8` | 110 / 10 | fail | pass | refined | 18.60 s | 3.61 s | 2.947 A |

Raw MPS passed 7/24. Refinement and validation-first selection passed 21/24:
11/12 gap-5 and 10/12 gap-10. Median selected backbone RMSD was 0.582 A and
median selected all-heavy RMSD was 1.352 A. Fixed-heavy RMSD and maximum
displacement were exactly zero for every published raw and refined candidate.

The three failures were `9ges`, `9g9t`, and `9e51`. The latter two also failed
the earlier Linux final-projection pilot. `9ges` passed that Linux baseline but
failed here with a generated/junction bond-length error after CPU refinement;
the same failure class and case had already appeared in the Linux per-step
reinjection follow-up. This is a device-sensitive case-level regression, not a
new failure class, and it remains in the denominator.

## Runtime And Memory

- Median MPS raw wall time: 19.14 s.
- Median CPU refinement wall time: 16.03 s.
- Maximum raw-plus-refinement wall time: 130.73 s (`24mc`).
- Linear-interpolation p95 raw-plus-refinement wall time: 92.47 s.
- Peak process RSS: 907,984,896 bytes (4.70% of physical unified memory).
- Peak sampled MPS tensor allocation: 233,721,344 bytes.
- Peak sampled Metal driver allocation: 355,254,272 bytes.

The memory values are overlapping unified-memory views and are not summed.
Raw MPS time remained close to 18–20 seconds from 19 to 461 residues, but the
CPU reference from the six-case gate was substantially faster for the shortest
case. Device-by-length performance policy remains unresolved.
