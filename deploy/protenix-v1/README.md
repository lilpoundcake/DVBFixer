# Protenix v1 Hook Spike

This directory contains the minimal maintained patch used to evaluate Protenix
v1 for reconstructing experimentally unresolved internal protein segments in
deposited structures. It is not a public DVBFixer backend.

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
