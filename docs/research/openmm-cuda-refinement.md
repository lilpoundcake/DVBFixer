# OpenMM CUDA Boundary-Refinement Environment

- Status: environment repaired; short-target repeatability and long-target gates pass.
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

OpenMM 8.2 predates the bundled `amber19/protein.ff19SB.xml` data used by
DVBFixer. The runtime and force-field data have independent compatibility
requirements: retain the OpenMM 8.2/CUDA 12.2 binary stack, then install the
unchanged XML files from the pinned OpenMM 8.3.0 source tag. The installer
verifies SHA-256 before modifying the environment:

```bash
bash deploy/openmm-cuda122/install-forcefields.sh /tmp/opencode/openmm-cuda122
```

The accepted digests are:

- `protein.ff19SB.xml`: `086bfdb1c05e7d5cb1d330e6c39fab5056a95e89492cdbd3044c746db13e9d09`;
- `tip3p.xml`: `3f4b188dbcb6c02863230eaca231e927fb6bf3307ce947d8a50d0f46f6dd83d9`.

Create it outside the repository and expose DVBFixer source explicitly:

```bash
micromamba create -p /tmp/opencode/openmm-cuda122 \
  -f deploy/openmm-cuda122/environment.yml
bash deploy/openmm-cuda122/install-forcefields.sh /tmp/opencode/openmm-cuda122
export PYTHONPATH="$PWD/src"
```

The local investigation environment currently lives at
`/tmp/opencode/openmm-cuda122`; that path is ephemeral and is not an artifact.

## Scientific Gates

The repaired environment successfully created an OpenMM CUDA Context on the A100
with `Precision=mixed` and `DeterministicForces=true`, evaluated energy, and
minimized a two-particle harmonic system. It reported OpenMM 8.2 and device
`NVIDIA A100-SXM4-40GB`, with final energy approximately zero.

The first invocation of the real adapter stopped during import because the
minimal environment initially omitted MDAnalysis. After adding it, the next
invocation exposed the absent OpenMM 8.2 ff19SB data described above. Neither
attempt reached a publishable scientific result.

With the repaired environment and exact ff19SB data, two independent `9dvi`
runs used explicit CUDA `Precision=mixed` and `DeterministicForces=true`. Both:

- passed every independent hard gate;
- preserved published fixed-heavy coordinates exactly (`0.0 Å` RMSD and maximum
  displacement);
- produced byte-identical PDB files with SHA-256
  `0284d877cbe4a5d4acd392aad549d5cbc09acac2172fdc53eede6cc7a3b6e945`;
- completed refinement in 3.31 and 1.79 seconds.

CUDA internally round-tripped massless fixed coordinates with a measured maximum
drift of `3.45e-6 Å`. Boundary-refinement revision v5 accepts only `<=1e-5 Å`
internal noise while publication continues to copy all fixed PDB records
unchanged. Missing identities or larger movement remain fatal.

The 440-residue `9ina` long-target gate then passed every hard gate in 13.83
seconds, with fixed-heavy RMSD/max displacement `0.0 Å`. Its published candidate
SHA-256 is `9a33bd3ff2b66ec5a4ba6483a943896edc0fa6bac2ce2b7289edc68448938aa1`.
The corresponding OpenMM Reference attempt had previously exceeded 900 seconds.

These gates authorize the full 231-case follow-up frozen in
[`small-diffusion-full-followup.json`](small-diffusion-full-followup.json),
SHA-256 `c9530bd95e1b420f27d31ae9c389940f096f66f215536f8626d2571689b6ba55`.
They do not establish cross-host or cross-GPU determinism.

The 231-case run will be descriptive rather than an independent confirmatory
benchmark because Protpardelle's exact training membership remains unavailable
and comparator outcomes are already known.
