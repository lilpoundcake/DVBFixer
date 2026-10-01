# Protpardelle-1c research adapter

This directory contains the pinned DVBFixer research adapter for Protpardelle-1c
revision `ee378400f25b801fa481028000f9060183d7fb4c`. It is not a production backend.

## Apple Silicon preparation

Use a native arm64 Python environment and apply the Apple patch to the clean,
pinned upstream checkout:

```bash
git -C "$PROTPARDELLE_ROOT" apply \
  "$DVBFIXER_ROOT/deploy/protpardelle-1c/apple-portability.patch"
```

The patch makes CUDA cache cleanup conditional and disables per-step host copies
of trajectory tensors when the adapter requests only final coordinates. It does
not change the denoising update or final coordinates. The existing callback patch
is not part of the frozen Apple baseline and must not be combined with this patch.

Before any MPS process starts, disable unsupported-operator fallback explicitly:

```bash
export PYTHONNOUSERSITE=1
export PYTORCH_ENABLE_MPS_FALLBACK=0
```

Freeze the remaining MPS variables rather than changing them between lanes. The
initial baseline leaves fast math and Metal preference disabled/unset and retains
PyTorch's default allocator watermarks:

```bash
unset PYTORCH_MPS_FAST_MATH
unset PYTORCH_MPS_PREFER_METAL
unset PYTORCH_MPS_HIGH_WATERMARK_RATIO
unset PYTORCH_MPS_LOW_WATERMARK_RATIO
```

Collect the host record from the Protpardelle environment before inference. The
repository must be a clean checkout and the environment export must already
exist:

```bash
python scripts/collect_apple_diffusion_inventory.py \
  .artifacts/apple-diffusion/host-inventory.json \
  --environment-lock .artifacts/apple-diffusion/protpardelle-explicit.txt
```

Before using a materialized DVBFixer request, run paired synthetic operator
smokes with the pinned config and checkpoint:

```bash
python deploy/protpardelle-1c/operator_smoke.py \
  --config "$CONFIG" --checkpoint "$CHECKPOINT" \
  --output .artifacts/apple-diffusion/cpu-operator-smoke.json --device cpu
python deploy/protpardelle-1c/operator_smoke.py \
  --config "$CONFIG" --checkpoint "$CHECKPOINT" \
  --output .artifacts/apple-diffusion/mps-operator-smoke.json --device mps
python deploy/protpardelle-1c/compare_operator_smokes.py \
  .artifacts/apple-diffusion/cpu-operator-smoke.json \
  .artifacts/apple-diffusion/mps-operator-smoke.json \
  .artifacts/apple-diffusion/cpu-mps-operator-comparison.json
```

Run the one-step MPS operator smoke with signpost profiling enabled:

```bash
python deploy/protpardelle-1c/checkpoint_gap_smoke.py \
  WORKSPACE/request.json WORKSPACE/apple-mps-one-step \
  --config "$CONFIG" --checkpoint "$CHECKPOINT" \
  --device mps --steps 1 --mps-profile
```

The adapter fails instead of relabeling CPU execution as MPS. It verifies all
model parameters and buffers, sampler inputs, and final coordinates are on the
requested device. Timings synchronize MPS before measurement boundaries. MPS
tensor allocation, Metal driver allocation, recommended working-set size, and
process peak RSS are recorded separately because they overlap in unified memory
and must not be added together.

OpenMM refinement runs in the separate DVBFixer environment with the `CPU`
platform. It is not MPS acceleration.
