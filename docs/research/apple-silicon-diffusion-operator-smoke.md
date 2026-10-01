# Apple Silicon Diffusion Operator Smoke

## Status

Preliminary portability evidence, collected before the frozen case-based Apple
pilot. This record proves that the pinned Protpardelle-1c checkpoint can execute
one complete denoising update on MPS without PyTorch CPU fallback. It does not
establish scientific utility, full-sampler latency, repeatability, or readiness
for integration.

## Frozen Inputs

- Protpardelle source revision:
  `ee378400f25b801fa481028000f9060183d7fb4c`
- Checkpoint: `cc89_epoch415.pth`
- Checkpoint SHA-256:
  `dfc9895b399ec4497bf6d646502168725f01dcbd55bc0b541fc3e3cc1f2f0483`
- Config SHA-256:
  `e9999ace79bf3044351cc982a624fc3459add3a4f91cbf97ff5b437b8942eb9d`
- Input: one synthetic 32-residue poly-alanine chain, seed `20260929`, one
  diffusion step, float32, supplied sequence, no MPNN.
- Apple patch: `deploy/protpardelle-1c/apple-portability.patch`. It guards an
  upstream CUDA-only cache call and disables unused per-step CPU trajectory
  copies. It does not change the denoising update or final coordinate tensor.

## Host And Environment

- Host: MacBook Pro `Mac15,6`, Apple M3 Pro, 11 CPU and 14 GPU cores.
- Unified memory: 19,327,352,832 bytes (reported as 18 GB by macOS).
- OS: macOS 26.6.2, arm64 native Python.
- Xcode Command Line Tools: `/Library/Developer/CommandLineTools`, Apple clang
  21.0.0 (`clang-2100.1.1.101`).
- Python: 3.12.14.
- PyTorch: 2.6.0, matching the Linux exploratory environment's framework
  version; MPS built and available.
- NumPy: 2.4.6.
- Environment installed size: 2,332,368 KiB.
- `PYTHONNOUSERSITE=1` and `PYTORCH_ENABLE_MPS_FALLBACK=0` were set before the
  smoke process started.
- `pip check`: no broken requirements.
- Ignored environment records:
  - `conda-explicit.txt` SHA-256
    `5ea8aab0f8c59f83168098a06df7b5a144ec0b508afc0a3c7c3172e65220e63d`;
  - `pip-freeze.txt` SHA-256
    `a83344873517dfc7a25aa652bffb0fbe6cacc9054a24717ff3a4f9c1956149b5`.

The upstream import warned that ESMFold, ProteinMPNN, LigandMPNN, and Foldseek
were absent. Those optional integrations are intentionally disabled by this
study and were not invoked.

## Result

| Measurement | CPU | MPS |
|---|---:|---:|
| Model load | 0.0711 s | 0.2158 s |
| One sampling step | 0.0314 s | 1.9739 s |
| Output | `(1, 32, 37, 3)`, finite | `(1, 32, 37, 3)`, finite on `mps:0` |
| MPS tensor allocation | n/a | 101,315,072 bytes |
| Metal driver allocation | n/a | 355,205,120 bytes |
| Recommended MPS working set | n/a | 14,302,248,960 bytes |

The paired output was not byte-identical. CPU-to-MPS coordinate RMSD was
`1.0341e-5 A`, mean displacement was `5.0055e-6 A`, and maximum displacement
was `1.2207e-4 A`. This is a small numerical backend difference, not a
scientific equivalence result.

A second independent MPS process with the same seed produced the same raw
coordinate SHA-256
`4ad30b255f6025c28708c44c985c9e8826c70431d214ce4dd3919cef54ac6e33`.
The two synthetic one-step MPS outputs are byte-identical. The repeat summary
SHA-256 is
`9e156890cbc9e7fa66a105c3bf5c9b6808565450c21d2bf0f1fd2372b4fbc5b7`.
This does not replace the planned independent-workspace 500-step `36hb`
repeatability test.

Ignored output records:

- CPU summary SHA-256:
  `683ba8d54e0d5970e3573adf9df350a4d45a9acf2cca3cf0b3a7d78d74c8b06d`;
- MPS summary SHA-256:
  `4832a0430e21b3b6a19516565dac34115fa80acb8f391e6025c19d435c98694b`;
- comparison SHA-256:
  `b537a01c9da62eecbe3ffa39bf34c8f2c4b425b998ce1809c72624ea17a050cc`.

The CPU result being faster for one short step does not predict 500-step
sampling performance. Cold Metal graph compilation, longer targets, repeated
steps, and thermal behavior remain unmeasured.

## Decision

The operator gate passes: the pinned all-atom checkpoint loads on native MPS,
the complete one-step sampler returns finite atom37 coordinates on MPS, and no
unsupported operator attempted PyTorch fallback. Proceed to the frozen `36hb`
500-step repeatability lane only after its workspace and comparator artifacts
are transferred and the clean-revision host inventory is frozen. Do not infer
case-level scientific suitability from this smoke.
