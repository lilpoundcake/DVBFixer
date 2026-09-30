# Diffusion Model Sizes And Evaluated Run Parameters

This record identifies the exact Protpardelle-1c, Protenix v1, and Boltz-2
models used in the DVBFixer gap-reconstruction studies. It separates model
parameter count, serialized checkpoint size, inference settings, and observed
runtime memory because those quantities are not interchangeable.

## Summary

| Engine | Evaluated model | Parameters | Checkpoint size | Evaluated schedule |
|---|---|---:|---:|---|
| Protpardelle-1c | `cc89_epoch415.pth` | 25,169,439 | 301,337,484 bytes (287.38 MiB) | 500 diffusion steps |
| Protenix v1 | `protenix_base_default_v1.0.0.pt` | 368.48M, upstream-reported | 1,475,950,125 bytes (1.375 GiB) | 1 cycle, 200 diffusion steps |
| Boltz-2 | `boltz2_conf.ckpt` | 506,724,992 | 2,286,561,469 bytes (2.130 GiB) | 1 recycling step, 200 sampling steps |

The Protpardelle model is approximately 14.6 times smaller than the evaluated
Protenix base model and 20.1 times smaller than the evaluated Boltz-2 model by
parameter count. Checkpoint-byte ratios differ because serialized files may
include metadata, duplicated/averaged weights, or other state; checkpoint size
must not be treated as parameter count or runtime memory.

## Protpardelle-1c `cc89`

- Source revision: `ee378400f25b801fa481028000f9060183d7fb4c`.
- Package version: 1.3.2.
- Checkpoint: `cc89_epoch415.pth`.
- Checkpoint SHA-256:
  `dfc9895b399ec4497bf6d646502168725f01dcbd55bc0b541fc3e3cc1f2f0483`.
- Checkpoint size: 301,337,484 bytes (301.34 MB; 287.38 MiB).
- Instantiated parameter count: 25,169,439.
- Config SHA-256:
  `e9999ace79bf3044351cc982a624fc3459add3a4f91cbf97ff5b437b8942eb9d`.

The evaluated `cc89` checkpoint is the all-atom model, not the backbone-only
`cc58` checkpoint. Its pinned config specifies a DiT structure model with an
atom37 output, 256 channels, 10 layers, two blocks per layer, eight attention
heads, a head dimension of 32, and a fixed maximum input size of 512 residues.

The frozen Apple Silicon 231-case run used:

- PyTorch 2.6.0 float32 MPS with `PYTORCH_ENABLE_MPS_FALLBACK=0`;
- 500 diffusion steps, `step_scale=1.2`, and `s_churn=200.0`;
- one candidate per request and the seed already frozen in that request;
- no MPNN sequence design and no per-step coordinate reinjection;
- final rigid synchronization followed by exact fixed-coordinate projection;
- generated-only OpenMM CPU refinement with eight restarts and no diffusion
  resampling.

Across that run, peak process RSS was 912,080,896 bytes, peak sampled MPS tensor
allocation was 289,383,168 bytes, and peak sampled Metal driver allocation was
355,254,272 bytes. These are overlapping unified-memory views and must not be
summed. See the
[full Apple results](apple-silicon-diffusion-full-231-results.md) for the
scientific and operational outcomes.

## Protenix v1 Base

- Source revision: `85767b811c40ed46e73a9b39519cf6bfca8701ba`.
- Checkpoint: `protenix_base_default_v1.0.0.pt`.
- Checkpoint SHA-256:
  `2b7d5a8b30494514fc47fd2271a16260528cdba170ba09cc112fdecd8f85ec04`.
- Checkpoint size: 1,475,950,125 bytes (1,475.95 MB; 1.375 GiB).
- Parameter count: 368.48M as reported by the
  [upstream supported-model record](https://github.com/bytedance/Protenix/blob/85767b811c40ed46e73a9b39519cf6bfca8701ba/docs/supported_models.md)
  for this exact model name.
- Training cutoff: 2021-09-30.

The confirmatory cohort used the base model, not Protenix Tiny or Mini. Each
case used one cycle, 200 diffusion steps, one candidate, template conditioning,
and exact fixed-atom reinjection after every denoising update. Generated-only
OpenMM refinement used the Reference platform and eight restarts without
rerunning diffusion.

On the Linux A100-SXM4-40GB confirmatory host, the complete per-case median was
209.0 seconds, peak host RAM was 3.68 GiB, and peak VRAM was 6.61 GiB. The
parameter count is an upstream model specification; the frozen DVBFixer record
independently verifies the checkpoint identity and byte size but did not publish
a second instantiated-parameter count.

## Boltz-2 Structure Model

- Source revision: `b1ebfc46ecf57f5414e0d1a6f9027bbb122c53bc`.
- Immutable Hugging Face revision:
  `6fdef46d763fee7fbb83ca5501ccceff43b85607`.
- Checkpoint: `boltz2_conf.ckpt`.
- Checkpoint SHA-256:
  `090e82ac8c92f5e943fa1b39e7410a44027bea7243c0bbb3caa67a77fc1428e1`.
- Checkpoint size: 2,286,561,469 bytes (2,286.56 MB; 2.130 GiB).
- Strictly loaded parameter count: 506,724,992.

The exact checkpoint byte size is also recorded by the
[pinned Hugging Face artifact](https://huggingface.co/boltz-community/boltz-2/blob/6fdef46d763fee7fbb83ca5501ccceff43b85607/boltz2_conf.ckpt).

Each confirmatory case used one recycling step, 200 sampling steps, one
diffusion sample, exact fixed-atom reinjection after every sampling update, and
one published candidate. Generated-only OpenMM refinement used the Reference
platform and eight restarts.

On the Linux A100-SXM4-40GB host, the complete per-case median was 160.0 seconds,
peak host RAM was 5.38 GiB, and peak VRAM was 3.93 GiB. The cache also contained
`boltz2_aff.ckpt` and `mols.tar`, but the gap-structure inference path loaded
`boltz2_conf.ckpt`; the affinity checkpoint was not the evaluated structure
model.

## Interpretation Limits

The three schedules are the frozen settings actually evaluated, not an
architecture-normalized speed benchmark. Protpardelle used 500 steps on Apple
MPS, whereas Protenix and Boltz-2 used 200-step schedules on an NVIDIA A100.
Their wall times and memory values therefore describe deployment behavior on
the tested hosts and must not be interpreted as model-size-only scaling.

For cohort outcomes and the eligibility distinction between Protenix and the
Boltz-2 proxy, see the
[confirmatory cohort results](diffusion-confirmatory-results.md). Artifact
origins, licenses, revisions, and distribution status remain recorded in
[`diffusion-gap-reconstruction-inventory.toml`](diffusion-gap-reconstruction-inventory.toml).
