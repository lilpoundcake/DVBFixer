# Protenix v1 Hook Spike

This directory contains the maintained patch, research drivers, and experimental
protocol runner used to evaluate Protenix v1 for reconstructing experimentally
unresolved internal protein segments in deposited structures. MODELLER remains
the default and production-supported backend.

Apply the patch only to revision
`85767b811c40ed46e73a9b39519cf6bfca8701ba`:

```bash
git -C /path/to/Protenix checkout --detach 85767b811c40ed46e73a9b39519cf6bfca8701ba
git -C /path/to/Protenix apply /path/to/DVBFixer/deploy/protenix-v1/per-step-callback.patch
python deploy/protenix-v1/hook-smoke.py /path/to/Protenix
```

The callback receives the mutable all-atom coordinate tensor after every
denoising update and before the next update. It must preserve tensor shape. The
benchmark adapter uses `reinjection.FixedAtomReinjector` and the stable atom
axis to perform weighted Kabsch synchronization on the active Torch device and
overwrite observed coordinates exactly.

The CPU smoke proves callback placement and exact overwrite mechanics without
loading weights. For the private checkpoint-backed gap smoke, build a template
input and run from an environment containing the pinned Protenix checkout:

```bash
python deploy/protenix-v1/build-template-input.py request.json protenix-input
PYTHONPATH=/path/to/Protenix:src \
  PROTENIX_ROOT_DIR=/path/to/Protenix \
  LAYERNORM_TYPE=torch \
  python deploy/protenix-v1/checkpoint_gap_smoke.py \
    request.json protenix-input/input.json protenix-output \
    --kalign /path/to/kalign \
    --checkpoint /path/to/protenix_base_default_v1.0.0.pt --steps 200
```

The driver resolves Protenix sequence ordinals back to request identities,
installs `FixedAtomReinjector` after featurization, records every callback, emits
only requested atoms, and runs DVBFixer-owned validation. The first full 7X35
run proved exact reinjection at all 200 steps but failed one peptide-junction
gate. Apply the existing DVBFixer localized refinement without rerunning the
checkpoint:

```bash
python deploy/protenix-v1/refine_candidate.py \
  request.json protenix-output/candidate.pdb refined-output \
  --expected-candidate-sha256 RAW_CANDIDATE_SHA256
```

The refinement wrapper preserves the raw candidate, rewrites only generated
coordinate columns, and reruns independent validation. Independent seeded 7X35
runs produced byte-identical raw and refined PDBs, and the refined candidate
passed every hard gate. The initial three-way ablation is complete: use
`--ablation-mode template-conditioning-only` for the conditioning-only arm;
the default is per-step reinjection. Broader benchmark evidence is still
required before any public backend work.

## Experimental CLI protocol runner

`production_runner.py` implements the versioned external-runner protocol for the
explicit `protenix-v1-cuda` profile. It freezes the audited Protenix revision,
callback patch, PyTorch/CUDA identity, 200-step per-step reinjection settings,
and deterministic OpenMM `Reference` boundary refinement.

The wrapper is committed and dependency-light tests run in core CI, but its
checkpoint-backed acceptance remains intentionally pending until it is executed
on a Linux/NVIDIA host. Do not treat the profile as hardware-accepted based on
macOS tests.

On the pinned Linux environment:

```bash
python "$DVBFIXER_ROOT/scripts/install_diffusion_backend.py" \
  protenix /path/to/protenix_base_default_v1.0.0.pt \
  --engine-root /path/to/Protenix

dvbfixer model INPUT.pdb --fasta TARGET.fasta \
  --backend diffusion -o OUTPUT_BUNDLE
```

The one-time installer verifies the frozen checkpoint digest and creates a
private launcher bound to the backend environment. Runtime profile, runner,
checkpoint, and digest are not user-facing model options.

The runner requires Linux amd64, Python 3.13, PyTorch 2.13.0, CUDA 12.9,
bfloat16-capable NVIDIA hardware, `kalign`, the exact patched source tree, and
the operator-provided checkpoint. DVBFixer does not download or redistribute
the checkpoint. Remove existing `__pycache__` directories from the Protenix
checkout before preflight; the wrapper disables bytecode writes before loading
the engine so accepted runs do not recreate them.
