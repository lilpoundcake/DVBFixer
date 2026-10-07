# Linux/NVIDIA Protenix Acceptance

- Evaluation date: 2026-10-07.
- Status: blocked before checkpoint loading; hardware acceptance is not claimed.
- Profile: `protenix-v1-cuda`.
- DVBFixer branch baseline: `bf6496de5e6da76db14282c2862847bb4406ef56`.

## Candidate Environment Identity

- Host GPU: NVIDIA A100-SXM4-40GB, compute capability 8.0, 40,960 MiB VRAM.
- Driver: `535.104.05`.
- Python: CPython `3.13.15`.
- PyTorch/CUDA: `2.13.0` / `12.9`.
- OpenMM: `8.6.1`.
- NumPy: `2.4.6`.
- Gemmi: `0.7.5`.
- PDBFixer: `1.12.0`.
- MDAnalysis: `2.10.0`.
- Protenix source revision:
  `85767b811c40ed46e73a9b39519cf6bfca8701ba`.
- Callback patch SHA-256:
  `cc4153be3dfd241124ea183d592884799300046b6ea7d3eccae8409b9fe21aa0`.
- Operator checkpoint SHA-256:
  `2b7d5a8b30494514fc47fd2271a16260528cdba170ba09cc112fdecd8f85ec04`.
- Kalign: `3.6.0`, executable SHA-256
  `057e7d91f5a56491e7dd097629b4b72a323295716a46de36d447581870597bbc`.
- Candidate profile-lock SHA-256:
  `b93eee004a4c10a1c7e1a55b2b9b25dcda914c361a7c03bfbee7efc82bdf621f`.
- Base image:
  `nvidia/cuda:12.9.1-cudnn-devel-ubuntu24.04@sha256:a2e1e2360c85298ac47ec2543b406ab1e8cec42e31ee47e4d32140ebc82e1067`.

The base-image manifest digest was independently resolved from the Docker Hub
registry. No OCI runtime was available on the host, so no final image was built
and no final image digest or complete package export is recorded. No checkpoint,
credential, private URL, or private path is tracked.

## Portable Runner Hardening

The production wrapper now fails closed unless all candidate profile identities
match before inference. It verifies the checkpoint before importing the adapter,
requires a bounded cgroup memory limit, rejects source changes outside the
maintained patch, verifies critical package and Kalign identities, checks the
exact A100 device class and CUDA runtime, and records fallback-disabled profile
facts.

The checkpoint adapter verifies model, sampler, and output device placement,
floating dtype continuity, contiguous callback indices and count, deterministic
algorithm state, finite output coordinates, and nonzero RAM/VRAM measurements.
The callback patch itself rejects shape, device, or dtype changes. Successful
runner results are accepted only when provenance, private artifact digests,
sampler trace, fallback policy, resource metrics, and refinement state all match
the frozen profile. Result JSON is written atomically and never replaces an
existing result.

## Preflight Result

The bounded preflight ran on the visible A100 without loading a model. It
reported:

- contract schema 4 and runner protocol 4;
- CUDA available on `cuda:0`;
- framework `2.13.0`;
- fallback disabled;
- expected patch and profile-lock identities;
- `missing-resource`: pinned Kalign was unavailable in the active environment;
- `missing-checkpoint`: the operator checkpoint was not provided;
- `incompatible-source`: the pinned Protenix checkout was unavailable.

This is the required fail-closed outcome. The operator-provided checkpoint path
was neither available nor recorded.

## Verification

- Focused Protenix, preflight, inventory, and sampler tests after final edits:
  `43 passed in 2.95s` with `PYTHONPATH=src`.
- Full non-slow repository suite with source isolation: `977 passed, 6 skipped,
  31 deselected in 260.59s`. The skips require unavailable GROMACS; existing
  Python 3.13 teardown warnings from `batch._Mirror` remain non-fatal.
- Ruff on the changed Python surfaces: passed.
- Callback patch application check against the exact source revision: passed.
- Base-image registry digest check: passed.
- `git diff --check`: passed.

## Stop Gate

Checkpoint-backed reconstruction, independent same-seed repeats,
different-seed candidates, failure injection, and the frozen 231-case replay
were not run. A self-hosted acceptance lane is not enabled until it can execute
the production wrapper for every outcome-bearing case rather than bypassing it
through a research adapter. No public scope, default, fallback, batch, ZBS, or
Homology claim changes. MODELLER remains the default.
