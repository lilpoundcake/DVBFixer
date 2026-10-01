# Boltz-2 Hook Spike

This directory contains the research-only maintained callback patch for Boltz-2
revision `b1ebfc46ecf57f5414e0d1a6f9027bbb122c53bc`. It is not a public DVBFixer
backend.

The patch exposes the mutable all-atom state immediately after each predictor
update and before the next denoising iteration. Apply and exercise it with:

```bash
git -C /path/to/boltz checkout --detach b1ebfc46ecf57f5414e0d1a6f9027bbb122c53bc
git -C /path/to/boltz apply /path/to/DVBFixer/deploy/boltz-2/per-step-callback.patch
PYTHONPATH=/path/to/boltz/src \
  python deploy/boltz-2/hook-smoke.py /path/to/boltz
```

The synthetic smoke invokes the real `diffusionv2.AtomDiffusion.sample` loop
for two steps and requires exact fixed-coordinate overwrite after each update.

The checkpoint-backed baseline uses the immutable Hugging Face revision
`6fdef46d763fee7fbb83ca5501ccceff43b85607`:

- `boltz2_conf.ckpt`: `090e82ac8c92f5e943fa1b39e7410a44027bea7243c0bbb3caa67a77fc1428e1`
- `boltz2_aff.ckpt`: `dcc5cd3722b1c9eaa34267e4ae32f55cbbf1963f4c19319381ccfa30fdd2ca9e`
- `mols.tar`: `39e076d96dbec6b4e86982bbda16f3a53a2a60c9bdc17828d88f6f9a0c7d1fd7`

Use a separate Python 3.12 environment. On the tested driver 535 host, pin a
CUDA 12.x PyTorch build rather than allowing current PyPI resolution to select
CUDA 13:

```bash
python -m pip install torch==2.5.1 \
  --index-url https://download.pytorch.org/whl/cu121
python -m pip install -e /path/to/boltz
python -m pip install numpy==1.26.4 openmm==8.6.1 MDAnalysis==2.10.0

BOLTZ_CACHE=/path/to/digest-verified-cache boltz predict \
  /path/to/boltz/examples/prot_no_msa.yaml \
  --model boltz2 --checkpoint /path/to/boltz2_conf.ckpt \
  --accelerator gpu --devices 1 --recycling_steps 1 --sampling_steps 2 \
  --diffusion_samples 1 --max_parallel_samples 1 --num_workers 0 \
  --no_kernels --seed 7 --output_format pdb --out_dir boltz-smoke
```

Two independent A100 runs produced byte-identical PDB and confidence outputs.
This stock command verifies the official model's real diffusion forward but
does not install a DVBFixer callback or provide constrained-gap evidence.

The research-only checkpoint gap path builds a complete-sequence query and a
deposited-coordinate template with explicit `SEQRES`, validates Boltz's token
and feature atom axes against exact DVBFixer identities, bypasses the stock
writer, and independently validates the candidate:

```bash
python deploy/boltz-2/build-template-input.py REQUEST.json INPUT_DIR

PYTHONPATH=/path/to/DVBFixer/src:/path/to/boltz/src \
  python deploy/boltz-2/checkpoint_gap_smoke.py \
  REQUEST.json INPUT_DIR/input.yaml OUTPUT_DIR \
  --cache /path/to/boltz-cache \
  --checkpoint /path/to/boltz2_conf.ckpt \
  --recycling-steps 1 --sampling-steps 200

python deploy/boltz-2/refine_candidate.py \
  REQUEST.json OUTPUT_DIR/candidate.pdb REFINED_DIR \
  --platform Reference --restart-count 8 \
  --expected-candidate-sha256 EXPECTED_SHA256
```

The 7X35 five-residue smoke produced byte-identical raw and refined repeats.
The raw candidate preserved fixed atoms exactly but missed one peptide-junction
gate; generated-only boundary refinement closed both junctions and passed every
frozen gate. Two diffusion steps exercise the callback but remain too immature
for candidate materialization: a generated coordinate exceeded PDB's fixed
column range. Use the 200-step schedule for scientific smoke evidence.
