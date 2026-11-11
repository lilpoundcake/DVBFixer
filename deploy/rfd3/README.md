# Private RFD3 Mac research hooks

This is **not a registered public DVBFixer backend** or a qualified reconstruction
profile. Execution evidence and limitations are recorded in
[`docs/research/mac-reconstruction-feasibility-2026-11-11.md`](../../docs/research/mac-reconstruction-feasibility-2026-11-11.md).

Reviewed source: `fnachon/foundry`
`7a9c246fbf140a743d14a35248f5039e7723a23a`. Its BSD-3-Clause license covers source
and referenced weights (Copyright 2025, Institute for Protein Design, University
of Washington); retain the upstream `LICENSE.md` with any redistribution. No
weights or assembled environment are shipped here.

- `scripts/trace_rfd3_seed_transport.py` records the **final model-consumer** fixed/
  sequence masks, RNG states, first noisy coordinates and native full-precision
  final array without changing the sampler. Optional expected-count environment
  variables fail before model forward on a mask mismatch.
- `scripts/rfd3_mps_seeded_research.py` defines the separately labeled
  `private-mps-seed-at-model-v1` arm. The untouched route's requested seeds reached
  different CPU states but identical MPS states/noisy arrays on this host. This
  entry explicitly seeds MPS at the model boundary, after loading/preprocessing.
  It does not establish the underlying reset cause or retroactively repair the
  original cohort.
- `source-component-width.patch` prevents source IDs such as `a1037` being
  truncated to `a103` by `set_indices()` using the chain-ID string dtype. Apply
  only to an owned clone. Record its digest and exact source diff. Earlier native
  context output remains intact; explicit consumer-bound metadata recovery is
  separate from resampling and does not prove all downstream feature equivalence.
- `scripts/prepare_rfd3_local_projection.py` uses the established request placement,
  zero missing-coordinate placeholders, complete canonical observations and
  exact chain case. Optional declared context spans use explicit `/0` cuts:
  these are synthetic model entities, **not original full-chain topology**.
  Crosslinks, unsupported insertions/numbering and unselected internal gaps fail
  closed. Other generated donors never enter conditioning.
- `scripts/materialize_rfd3_local_controls.py` checks native sequence/inventory,
  explicit index maps, observed-only frame fitting and generated ownership before
  inserting a group into a separately retained full assembly. It does not invent
  a Protpardelle `SamplerTrace` or grant publication rights.

Synthetic isolated-runtime tests (not checkpoint evidence):

```bash
PYTHONPATH="$PWD/scripts:$PWD/src" PYTORCH_ENABLE_MPS_FALLBACK=0 \
  /path/to/owned/rfd3/bin/python -m pytest deploy/rfd3/test_mps_research_hooks.py -q
```

All real runs use the serial ML lock and `pressure-v2` guard in
`scripts/run_mac_reconstruction_guarded.py`. Do not disable MPS watermarks,
implicitly crop a failed request or substitute arbitrary checkpoints. Source
preservation, raw quality, derivative refinement and broad qualification remain
distinct gates.
