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

Unmodified Boltz-2 remains unsupported. DVBFixer now carries the minimal
Boltz-2-only patch at `deploy/boltz-2/per-step-callback.patch`. It exposes the
mutable state after `atom_coords = atom_coords_next`, validates the callback's
type, shape, device, and dtype, and threads the callback through `Boltz2`.
`deploy/boltz-2/hook-smoke.py` invoked the real sampler loop for two updates on
both CPU and A100; both callbacks restored fixed coordinates with `0.0 Å`
error and the returned final state retained them exactly.

The immutable official structure checkpoint has SHA-256
`090e82ac8c92f5e943fa1b39e7410a44027bea7243c0bbb3caa67a77fc1428e1`.
It loaded strictly as a 506,724,992-parameter `Boltz2` model. Two independent
seed-7 A100 runs of the upstream 120-residue, empty-MSA example with one recycle
and two diffusion steps produced byte-identical PDBs (`1ef1df41...b61b7`) and
zero failed examples. This stock-CLI baseline did not exercise the callback.
Stable Boltz atom-axis to DVBFixer identity mapping and checkpoint-backed gap
reinjection remain unresolved, so Boltz is not a constrained-gap backend.

## Checkpoint-Backed Gap Smoke

`deploy/protenix-v1/checkpoint_gap_smoke.py` exercised the official v1
checkpoint on the withheld five-residue `7x35-chain-a-glypro-5` case. Protenix
featurized 267 residues into a stable 2,118-atom axis; the request mapping
resolved all 2,117 fixed/generated output atoms without ambiguity. On an A100,
seed 7 with one cycle and 200 diffusion steps completed in 96.46 seconds. Every
post-update callback reported exactly `0.0 Å` fixed-coordinate error. Seeding
before dataloader construction is required because featurization creates random
reference conformers. With that ordering, two independent runs produced the
same `e1c7ecf...aa3e` PDB byte-for-byte.

Independent DVBFixer validation found `0.0 Å` fixed-heavy RMSD and maximum
displacement, complete generated residues and canonical heavy atoms, no
detectable D-Cα centers, no severe overlaps, and no local bond/angle/amide/
Ramachandran/chi failures. The raw candidate did not pass overall: its two
junction C-N distances were `1.284 Å` and `1.600 Å`, so the second exceeds the
frozen `1.45 Å` maximum. This is positive hook evidence and negative raw-output
evidence; it is not a publishable reconstruction result.

`deploy/protenix-v1/refine_candidate.py` then applied the existing pinned
`dvbfixer-openmm-boundary-refinement-v3` policy to generated atoms only. The
refined candidate passed every hard gate with junctions at `1.339 Å` and
`1.348 Å`, fixed-heavy RMSD `0.0 Å`, no D-Cα centers or severe overlaps, and
gap-backbone RMSD `1.964 Å`. Refining both independent checkpoint outputs
produced the same `38f0195c...490f` PDB byte-for-byte.

## Three-Way Ablation

The pinned checkpoint, case, seed, cycle count, and 200-step schedule were held
constant. Template conditioning alone used only one final Kabsch frame
synchronization: both peptide junctions passed, but fixed-heavy RMSD was
`1.986 Å` and maximum displacement was `10.263 Å`. Adding exact per-step
reinjection reduced both fixed-coordinate metrics to `0.0 Å`, but the right
junction remained open at `1.600 Å`. Adding generated-only localized boundary
refinement retained `0.0 Å` fixed-coordinate movement and closed the junctions
to `1.339 Å` and `1.348 Å`; every frozen hard gate passed.

This completes the predeclared three-way ablation for the initial 7X35 case.
It demonstrates that template conditioning is insufficient, per-step
reinjection is necessary for exact framework preservation, and localized
refinement is necessary for peptide closure on this sample. It is not broad
corpus or production-backend acceptance evidence.

## Additional Protenix Cases

The same seed-7, 200-step reinjection and generated-only refinement path was
repeated in independent workspaces for three additional target-chain-only
cases. Raw and refined PDBs were byte-identical between repeats in every case,
fixed-heavy RMSD and maximum displacement remained `0.0 Å`, and all refined
candidates passed every frozen hard gate:

| Case | Gap | Raw result | Refined junctions | Backbone RMSD |
|---|---:|---|---|---:|
| 8CZ8 withheld loop | 10 | passed | `1.329/1.340 Å` | `0.175 Å` |
| 7X35 Pro/Gly-rich | 7 | failed short `1.140 Å` junction | `1.337/1.339 Å` | `0.562 Å` |
| 7K8S H/82A-H/82C | 3 | passed | `1.330/1.344 Å` | `1.351 Å` |

The insertion-code identities survive end to end. The 7X35 result independently
confirms that localized refinement is an active closure step rather than a
no-op. Multichain interface and retained-heterogen cases remain outside this
single-chain Protenix driver.

## Decision

Unmodified Protenix and Boltz remain unsupported because they expose no required
hook. The maintained Protenix v1 patch plus DVBFixer localized boundary
refinement passes the initial-case three-way ablation and three additional
repeatable single-chain cases. The maintained Boltz-2 patch passes CPU/GPU hook
smokes and its official checkpoint passes repeatable real diffusion inference,
but Boltz request-identity mapping is still absent. Multichain/chemical-context,
environment, and distribution work remain open. Public CLI work must not begin
from these results.

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
