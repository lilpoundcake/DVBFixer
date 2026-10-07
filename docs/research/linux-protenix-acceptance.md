# Linux/NVIDIA Protenix Acceptance

- Evaluation date: 2026-10-07.
- Status: single-case hardware acceptance passed; full WP1 acceptance remains
  blocked on the frozen 231-case replay and self-hosted lane.
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
  `3fb3661be2b6748650b90bf9b0aaa7a09cbba9accc13253a6c29f644001b5d2e`.
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
maintained patch, verifies exact patched-source and auxiliary-artifact digests,
verifies critical package and Kalign identities, checks the exact A100 device
class and CUDA runtime, and records fallback-disabled profile facts.

The checkpoint adapter verifies model, sampler, and output device placement,
floating dtype continuity, contiguous callback indices and count, deterministic
algorithm state, finite output coordinates, and nonzero RAM/VRAM measurements.
The callback patch itself rejects shape, device, or dtype changes. Successful
runner results are accepted only when provenance, private artifact digests,
sampler trace, fallback policy, resource metrics, and refinement state all match
the frozen profile. Result JSON is written atomically and never replaces an
existing result.

## Preflight Result

The bounded preflight ran on the visible A100 with the exact local checkpoint,
source tree, auxiliary artifacts, and Kalign executable. It reported:

- contract schema 4 and runner protocol 4;
- CUDA available on `cuda:0`;
- framework `2.13.0`;
- fallback disabled;
- expected patch and profile-lock identities;
- exact checkpoint, patched-source, and auxiliary-artifact identities;
- no issues.

The tracked evidence records only the checkpoint digest, never its private path.
An injected invalid checkpoint used `README.md` as the input and failed closed
with only `missing-checkpoint` and `checkpoint-digest-mismatch`; no model was
loaded and no result was published.

## Checkpoint-Backed Reconstruction

Case `36hb` was rematerialized through the current schema-4 confirmatory cohort
builder and executed twice in independent workspaces with seed 7. Both refined
candidates passed every independent hard gate, preserved fixed heavy atoms at
`0.0` Angstrom RMSD and maximum displacement, recorded 200 callback updates with
maximum post-projection error `0.0` Angstrom, and produced the byte-identical
candidate SHA-256
`664e65f90a9a7ec6abdd762008b6f28fe20d2e47eaab576bdf1830a246b61ea4`.
Sampler traces were identical after excluding per-run RAM and wall-time fields.

The two seed-7 runs used 121.39 and 118.59 seconds, respectively, peaked at
3,051,652,096 bytes VRAM, and peaked at 3,935,883,264 and 3,915,182,080 bytes
RAM. A separately materialized seed-11 run passed every hard gate and produced
the distinct candidate SHA-256
`effcb40bbd42fdba2f95924bd51d3e851cac6f4f594e64f9240ae8d303dd1af4`.
It used 159.98 seconds, 3,051,652,096 bytes peak VRAM, and 3,926,147,072 bytes
peak RAM.

The production smoke also exposed and fixed a publication regression: Protenix's
canonical tensor may include model-only terminal `OXT`, while the candidate and
sampler-trace contracts permit exactly the represented fixed and requested
generated atoms. The adapter now projects the canonical axis onto that exact
request set before materialization and trace publication, without weakening the
full-tensor per-step reinjection callback.

## Public CLI Acceptance

The installed private launcher was exercised through the public command, not by
calling the research adapter directly:

```bash
dvbfixer model normalized.pdb --fasta target.fasta \
  --backend diffusion --diffusion-model protenix \
  --diffusion-work-parent /tmp/opencode -o /tmp/opencode/protenix-cli-e2e-36hb
```

Case `36hb` completed request construction, scope admission, runner preflight,
subprocess execution, A100 sampling, localized refinement, independent
validation, and atomic bundle publication. The bundle status is `success`, its
candidate has no hard-gate failures, and the final PDB SHA-256 is the same
deterministic seed-7 identity recorded above:
`664e65f90a9a7ec6abdd762008b6f28fe20d2e47eaab576bdf1830a246b61ea4`.

This run found and fixed two integration defects that direct adapter execution
could not expose: the installed launcher now binds the current DVBFixer source,
pinned engine source, backend Python, and digest-verified Kalign; core preflight
now expects the production runner's exact profile-lock environment identity.
The repository bind mount does not support `renameat2(RENAME_NOREPLACE)`, so the
accepted run used a native `/tmp/opencode` work/output filesystem. Publication
remains fail-closed on filesystems without atomic no-replace directory rename.

## Verification

- Focused installer, Protenix preflight, doctor, provenance, CLI, runtime, and
  pipeline tests after final edits: `69 passed in 2.00s` with `PYTHONPATH=src`.
- Full non-slow repository suite with source isolation: `990 passed, 6 skipped,
  31 deselected in 259.90s`. The skips require unavailable GROMACS; existing
  Python 3.13 teardown warnings from `batch._Mirror` remain non-fatal.
- Ruff on the changed Python surfaces: passed.
- Callback patch application check against the exact source revision: passed.
- Base-image registry digest check: passed.
- `git diff --check`: passed.

## Stop Gate

Checkpoint-backed reconstruction, independent same-seed repeats, a
different-seed candidate, and checkpoint failure injection passed. The frozen
231-case production-wrapper replay has not run, no final OCI image digest is
available, and a self-hosted acceptance lane is not enabled. Therefore the
accepted baseline has not yet been reproduced end to end and full WP1 acceptance
is not claimed. No public scope, default, fallback, batch, ZBS, or Homology claim
changes. MODELLER remains the default.
