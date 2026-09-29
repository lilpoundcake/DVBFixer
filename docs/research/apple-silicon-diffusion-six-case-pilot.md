# Apple Silicon Protpardelle Six-Case Pilot

## Status

The frozen six-case Apple Silicon gate passed on 2026-09-30. Protpardelle-1c
completed every 500-step sample natively on MPS with PyTorch CPU fallback
disabled. Generated-only OpenMM CPU refinement produced six of six candidates
that passed every DVBFixer hard gate. This is feasibility evidence for compact
diffusion inference on Apple Silicon, not evidence that MPS is faster than the
native CPU.

## Frozen Setup

- Host: Apple M3 Pro, 14 GPU cores, 19,327,352,832 bytes unified memory,
  macOS 26.6.2, arm64.
- Protpardelle source revision:
  `ee378400f25b801fa481028000f9060183d7fb4c`.
- PyTorch 2.6.0, float32, `PYTORCH_ENABLE_MPS_FALLBACK=0`.
- Checkpoint SHA-256:
  `dfc9895b399ec4497bf6d646502168725f01dcbd55bc0b541fc3e3cc1f2f0483`.
- Config SHA-256:
  `e9999ace79bf3044351cc982a624fc3459add3a4f91cbf97ff5b437b8942eb9d`.
- Sampling: 500 steps, step scale 1.2, churn 200.0, request seed 7,
  supplied sequence, MPNN disabled, one candidate, no per-step reinjection.
- Cases: the preregistered `36hb`, `9gtp`, `9gae`, `9ina`, `9eho`, and
  `9gbg` set. Workspaces were deterministically rematerialized from the tracked
  frozen cohort metadata and exact RCSB mmCIF URLs.
- Refinement: generated-only OpenMM CPU, eight restarts, fixed atoms projected
  exactly. Refinement is not MPS acceleration.

The partial materialization report marks cases as `awaiting-full-screen` because
the other 494 cohort cases were intentionally not downloaded. Membership in
this gate comes from the preregistered six-case Apple plan, not from changing
selection after observing results.

## Portability Finding

The request-level one-step smoke exposed one post-denoising portability defect:
the adapter attempted to widen an MPS tensor to float64 before copying it to
CPU. MPS does not implement float64. The adapter now performs
`MPS float32 -> CPU -> float64 NumPy`; a dependency-light regression test locks
that order. The denoiser itself had already completed on MPS, and no fallback
was enabled or used.

## Repeatability And CPU Reference

Two independent 500-step MPS processes for `36hb` passed all hard gates and
produced byte-identical candidate PDBs:

- candidate SHA-256:
  `adf0acfd6bc010eae01c8c493a1a698d0ae1f8a6d19fed76c9dcf7e8e7405528`;
- generated-atom RMSD, mean displacement, and maximum displacement: all
  `0.0 A`;
- wall times: 18.149 s and 18.336 s;
- peak RSS: 892,436,480 and 894,959,616 bytes.

The full 500-step CPU reference also passed. It was faster on this short target:
5.890 s total (5.807 s denoising) versus 18.149 s total (17.864 s denoising)
on MPS, a 3.08x MPS slowdown. CPU and MPS outputs were both valid but were not
expected to match: generated-atom RMSD was 1.539 A and maximum displacement was
4.528 A. Their gap-backbone RMSDs to the withheld experimental coordinates were
0.226 A (CPU) and 0.245 A (MPS).

## Six-Case Results

| Case | Length / gap | Raw MPS | Raw failure | MPS wall | Raw backbone RMSD | Refined CPU | Refine wall | Refined backbone RMSD |
|---|---:|---:|---|---:|---:|---:|---:|---:|
| `36hb` | 106 / 5 | pass | — | 18.149 s | 0.245 A | pass | 4.566 s | 0.241 A |
| `9gtp` | 225 / 5 | pass | — | 19.367 s | 0.467 A | pass | 11.767 s | 0.499 A |
| `9gae` | 103 / 5 | pass | — | 18.416 s | 0.496 A | pass | 6.204 s | 0.342 A |
| `9ina` | 440 / 10 | fail | severe overlap | 19.517 s | 0.522 A | pass | 35.806 s | 0.378 A |
| `9eho` | 118 / 10 | fail | severe overlap | 18.641 s | 1.046 A | pass | 4.241 s | 0.981 A |
| `9gbg` | 287 / 10 | fail | junction connectivity | 19.628 s | 1.517 A | pass | 21.286 s | 1.389 A |

Raw MPS passed 3/6. Validation-first selection after the frozen refinement path
passed 6/6. Fixed-heavy-atom RMSD and maximum displacement were exactly 0.0 A
for every published raw and refined candidate. Median MPS wall time was 19.004 s;
median raw backbone RMSD was 0.509 A. Median refinement wall time was 8.986 s;
median refined backbone RMSD was 0.438 A.

Peak observed resource values across the six independent MPS processes were:

- process RSS: 909,770,752 bytes (4.71% of physical unified memory);
- sampled MPS tensor allocation: 231,621,632 bytes;
- sampled Metal driver allocation: 355,336,192 bytes.

These are overlapping unified-memory views and are not summed.

## Candidate Digests

| Case | Raw MPS candidate SHA-256 | Refined candidate SHA-256 |
|---|---|---|
| `36hb` | `adf0acfd6bc010eae01c8c493a1a698d0ae1f8a6d19fed76c9dcf7e8e7405528` | `792d02c8732c055db0f11587801846d6e8bb58583c23f5def8fdc496b18bdbb0` |
| `9gtp` | `d47f4cd94bf229454adf327f6995b22c951b2b33947732bfb35570ee9092a176` | `5b522335a4e27134c209a61972b4ad589536e588a51d3a4d1fcbbf75388570b1` |
| `9gae` | `1fdfefa42acebf2d4e65339614e742af7f087d8773523af890e00a4c937fe881` | `69342c2bd44f29bc614ae56c43211ec390e27e4b5314059abbe191cea6bd44fa` |
| `9ina` | `72f757ac835d0cf3a084e2d0385b99aaf924f5b6d3964040467b517e8e4e0dc7` | `c7cbf0d700c7114c36d0e756878fdd6301c32099c1aca54729a7510e44d6480e` |
| `9eho` | `a27371e08b1c27a7f9afd1be3a1db13187ddd9a90996635de8c85b59701407c5` | `77ffc81945159844e16b76045b3edfab6813fd12b59bf485f6ce195c5a761949` |
| `9gbg` | `708989c065c2a71ca475a4d6f1058d1823791a4837dca519e2ab570a41acb78c` | `1beec2e59c727eefd54897ab39b17aad13f8f6fa958545540450e0417a4fc98b` |

## Decision

The six-case gate permits proceeding to the frozen 24-case portability pilot.
MPS is operationally viable and comfortably within memory limits, but the
current M3 Pro result gives no performance reason to prefer MPS over CPU for the
short `36hb` workload. The 24-case lane must retain both device/runtime data and
scientific failures in its denominator before deciding whether MPS, CPU, or a
device-by-length policy is appropriate on Apple Silicon.
