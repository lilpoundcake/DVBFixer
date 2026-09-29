# OpenMM CUDA Boundary-Refinement Environment

- Status: environment repaired; scientific CUDA refinement tests pending.
- Diagnosis date: 2026-09-29.
- Host GPU: NVIDIA A100-SXM4-40GB.
- NVIDIA driver: `535.104.05`.
- Driver-advertised CUDA compatibility: `12.2`.

## Root Cause

The main `dvbfixer` environment resolved OpenMM 8.6.1 together with CUDA/NVRTC
12.9 (`cuda-nvrtc 12.9.86`, `cuda-version 12.9`). OpenMM's CUDA platform JIT
therefore emitted PTX from the CUDA 12.9 toolchain. The installed 535-series
driver supports CUDA 12.2 and rejected that PTX before refinement with:

```text
CUDA_ERROR_UNSUPPORTED_PTX_VERSION (222)
```

This was not a PyTorch or A100 failure. Protpardelle sampling used its isolated
PyTorch CUDA 12.4 environment successfully. The failure was specific to the
OpenMM plugin loading `libnvrtc.so.12` from the main environment's CUDA 12.9
runtime.

## Isolated Repair

The reproducible environment specification is
[`../../deploy/openmm-cuda122/environment.yml`](../../deploy/openmm-cuda122/environment.yml).
It deliberately leaves the main scientific environment unchanged and pins:

- Python 3.12.14;
- OpenMM 8.2.0;
- CUDA version 12.2;
- NVRTC 12.2.140;
- cuFFT 11.0.8.103;
- PDBFixer 1.12;
- the dependencies needed by DVBFixer validation.

Create it outside the repository and expose DVBFixer source explicitly:

```bash
micromamba create -p /tmp/opencode/openmm-cuda122 \
  -f deploy/openmm-cuda122/environment.yml
export PYTHONPATH="$PWD/src"
```

The local investigation environment currently lives at
`/tmp/opencode/openmm-cuda122`; that path is ephemeral and is not an artifact.

## Evidence Before Scientific Tests

The repaired environment successfully created an OpenMM CUDA Context on the A100
with `Precision=mixed` and `DeterministicForces=true`, evaluated energy, and
minimized a two-particle harmonic system. It reported OpenMM 8.2 and device
`NVIDIA A100-SXM4-40GB`, with final energy approximately zero.

The first invocation of the real DVBFixer refinement adapter stopped during
Python import because the minimal environment initially omitted MDAnalysis. It
did not reach OpenMM system creation and is not a scientific CUDA result.
MDAnalysis 2.10.0 and Biopython 1.88 were subsequently added and are included in
the tracked environment specification.

No real candidate has been refined with the repaired stack at this documentation
checkpoint. The next required gates, in order, are:

1. Run one fixed `9dvi` refinement twice on OpenMM CUDA.
2. Require both outputs to pass the existing independent geometry validation.
3. Measure coordinate repeatability rather than inferring it from
   `DeterministicForces=true`.
4. Run one longer target to expose memory or kernel-compilation regressions.
5. Freeze a full 231-case follow-up manifest and run only after these gates pass.

The 231-case run will be descriptive rather than an independent confirmatory
benchmark because Protpardelle's exact training membership remains unavailable
and comparator outcomes are already known.
