# Linux Diffusion Expansion Baseline

- Inventory date: 2026-10-07.
- Branch tip: `bf6496de5e6da76db14282c2862847bb4406ef56`.
- Scope: WP0 baseline only; this is not Linux Protenix acceptance evidence.
- Source version: `0.9.0`.

## Host Inventory

- OS/kernel: Linux `5.15.0-102-generic`, `x86_64`.
- CPU: 2 x AMD EPYC 7742, 128 available physical cores.
- RAM: 1.5 TiB available to the host; no swap.
- Runner cgroup memory limit: 138,512,695,296 bytes.
- Local workspace storage allocation: 512 GiB total, 396 GiB available at
  inventory time. The private mount path is intentionally not recorded.
- GPU: NVIDIA A100-SXM4-40GB, compute capability 8.0, 40,960 MiB VRAM.
- Driver: `535.104.05`; driver-reported CUDA compatibility: `12.2`.
- Active Torch runtime: `2.13.0`, CUDA `12.9`; an allocation and readback on
  `cuda:0` succeeded and reported the expected A100 identity.
- Python: CPython `3.13.15`.
- Scientific packages: NumPy `2.4.6`, OpenMM `8.6.1`, PDBFixer `1.12.0`,
  Gemmi `0.7.5`, MDAnalysis `2.10.0`.
- OCI/runtime tools: Docker, Podman, Singularity, Apptainer, and
  `nvidia-container-cli` were not available.
- Protenix source and `kalign`: not available in the active environment.
- Operator checkpoint: not provided. No checkpoint path is recorded.

The active installed `dvbfixer` distribution metadata reported `0.8.9`, while
the checked-out source and synchronized release metadata report `0.9.0`.
The initial baseline test invocation did not set `PYTHONPATH=src`; a later
warning proved that this editable installation could resolve modules from an
older worktree. The recorded initial count below is therefore inventory only,
not authoritative validation of the frozen checkout. Final WP1 verification
was repeated with explicit source isolation. A production runner environment
must install the exact accepted source rather than reuse this developer
environment implicitly.

## Frozen Scientific Baseline

- Diffusion contract schema: 4.
- External runner protocol: 4.
- Protenix profile status: `hardware-acceptance-pending`.
- Protenix source revision:
  `85767b811c40ed46e73a9b39519cf6bfca8701ba`.
- Checkpoint SHA-256:
  `2b7d5a8b30494514fc47fd2271a16260528cdba170ba09cc112fdecd8f85ec04`.
- Fallback policy: disabled across profile, device, and MODELLER.
- Fixed-heavy RMSD maximum: 0.01 A.
- Fixed-heavy maximum displacement: 0.03 A.
- Junction C-N range: 1.20-1.45 A.
- Generated-backbone break maximum: 1.80 A.
- Severe clash overlap threshold: 0.90 A.

The frozen confirmatory evidence remains the 231-case v5 cohort documented in
`diffusion-confirmatory-results.md`: Protenix minimized passed 220/231
(95.24%, Wilson 95% CI 91.68-97.32%) versus MODELLER 169/231 (73.16%). The
tracked aggregate identity is
`55e4cbc6c16ee6bf34a14aabf61ec71afcddfb436f15edb5195c3d8b3f717ad7`.
The ignored aggregate itself was not present in this clean worktree, so this
inventory freezes the existing tracked identity and does not claim a new
aggregation run.

## Baseline Verification

Commands were run from a second clean worktree at the exact branch tip before
WP1 changes:

- Initial CPU diffusion suite: `258 passed in 236.82s`, but non-authoritative
  because of the editable-install contamination described above.
- `ruff check src/dvbfixer`: passed.
- `python scripts/check_versions.py`: passed at `0.9.0`.
- `python scripts/gen_gui_spec.py --check`: passed.
- `python scripts/check_agent_docs.py`: passed.
- `python scripts/gen_cli_reference.py --check`: failed with 22 drifted command
  pages under Python 3.13. Regeneration showed only the known Python 3.13
  argparse option-invocation formatting change; repository policy records that
  this check is generated under Python 3.11, which was unavailable on this
  host. No generated references were changed.
- Strict mypy command from `AGENTS.md`: failed with 25 errors in 10 files. The
  errors are pre-existing at the frozen tip and include diagnose geometry and
  chemistry typing, batch stream types, biological-assembly tuple widths,
  SMILES adjacency annotations, diffusion heterogen-context tuple widths,
  antibody list inference, and model pipeline option types.

WP1 checkpoint inference, repeatability, failure injection, and the 231-case
production-wrapper replay remain blocked until the operator-provided checkpoint,
pinned Protenix source, `kalign`, and a complete accepted environment are
available. This baseline does not widen any public scientific scope.
