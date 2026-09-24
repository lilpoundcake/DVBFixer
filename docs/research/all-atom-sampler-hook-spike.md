# All-Atom Sampler Hook Spike

This note records the Phase 3 source audit and checkpoint-backed hook smoke
performed on 2026-09-24. It establishes constrained-sampling control, not
production model quality or public-backend readiness.

## Required Hook

The adapter must receive the mutable all-atom coordinate state after every
denoising update and before the next update, retain an exact mapping from the
tensor atom axis to DVBFixer atom identities, and overwrite every fixed atom
exactly. Ordinary template conditioning and energy guidance are insufficient.

## Protenix v1

Audited revision: `85767b811c40ed46e73a9b39519cf6bfca8701ba`.

`protenix/model/generator.py:sample_diffusion` owns `x_l` inside its sampling
loop and returns it only after all steps. It accepts no state callback. The TFG
path calls `protenix/tfg/engine.py:TFGEngine.step`, but its projection modifies
the denoised `x0` estimate before the predictor-corrector update. The returned
`x_next` is not subsequently passed through an externally supplied exact
fixed-coordinate overwrite.

The unmodified source is unsupported. Monkeypatching `TFGEngine.step` is not a
stable upstream API. DVBFixer therefore carries a minimal patch at
`deploy/protenix-v1/per-step-callback.patch` that adds an optional callback
after each `x_l` update and before the next iteration. The patch applies cleanly
to the pinned revision. `deploy/protenix-v1/hook-smoke.py` exercised the real
patched `sample_diffusion` function for two synthetic CPU steps and observed
two exact fixed-coordinate overwrites. The callback implementation in
`deploy/protenix-v1/reinjection.py` performs weighted Kabsch synchronization on
the active Torch device, validates the stable `AtomIdentity` axis, and then
overwrites every fixed coordinate exactly.

## Boltz-2

Audited revision: `b1ebfc46ecf57f5414e0d1a6f9027bbb122c53bc`.

`src/boltz/model/modules/diffusionv2.py:AtomDiffusion.sample` owns
`atom_coords` inside its loop and exposes no state callback. Physical/contact
steering changes `atom_coords_denoised`; the sampler then computes
`atom_coords_next` and assigns it to `atom_coords` without an external
post-update hook. The legacy loop in
`src/boltz/model/modules/diffusion.py:AtomDiffusion.sample` has the same shape.

Result: unsupported. Boltz-2 also requires an upstream callback or maintained
source patch before exact per-step reinjection can be demonstrated.

## Checkpoint-Backed Gap Smoke

`deploy/protenix-v1/checkpoint_gap_smoke.py` exercised the official v1
checkpoint on the withheld five-residue `7x35-chain-a-glypro-5` case. Protenix
featurized 267 residues into a stable 2,118-atom axis; the request mapping
resolved all 2,117 fixed/generated output atoms without ambiguity. On an A100,
seed 7 with one cycle and 200 diffusion steps completed in 96.19 seconds. Every
post-update callback reported exactly `0.0 Å` fixed-coordinate error.

Independent DVBFixer validation found `0.0 Å` fixed-heavy RMSD and maximum
displacement, complete generated residues and canonical heavy atoms, no
detectable D-Cα centers, no severe overlaps, and no local bond/angle/amide/
Ramachandran/chi failures. The raw candidate did not pass overall: its two
junction C-N distances were `1.254 Å` and `1.604 Å`, so the second exceeds the
frozen `1.45 Å` maximum. This is positive hook evidence and negative raw-output
evidence; it is not a publishable reconstruction result.

## Decision

Unmodified Protenix and Boltz remain unsupported because they expose no required
hook. The maintained Protenix v1 patch is selected for the next localized
boundary-refinement spike because checkpoint-backed identity mapping and exact
per-step reinjection now pass. Boundary refinement, the three-way ablation,
repeatability, broader benchmarks, and environment/distribution work remain
open. Public CLI work must not begin from this result.

The pinned source publishes the v1 checkpoint URL as
`https://protenix.tos-cn-beijing.volces.com/checkpoint/protenix_base_default_v1.0.0.pt`.
The host's configured HTTPS proxy failed against that endpoint, while a direct
connection to the same official URL succeeded. The complete 1,475,950,125-byte
file has independently measured SHA-256
`2b7d5a8b30494514fc47fd2271a16260528cdba170ba09cc112fdecd8f85ec04`.
A callback-backed gap-reconstruction run has now completed; no mirror was used.

A separate base-model smoke loaded that checkpoint strictly on an A100 40 GB
and completed one 20-residue, 168-atom sample with one cycle and two diffusion
steps. Model forward time was 3.14 seconds. This confirms checkpoint/runtime
compatibility only: the upstream CLI did not pass the maintained callback, so
the run provides no gap-reconstruction or fixed-coordinate evidence.
