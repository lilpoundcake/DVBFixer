# Private fixed-source independent-cohort validation

Date: 2026-11-11. User authorized autonomous continuation to completion or a real
decision boundary. This freezes a **zero-observed-motion** research pipeline,
not the separate positive anchor-mobility experiment or a public profile.

## Scope and independent cases

Use the prospectively named 1PGA (protein-G B1, A1–56) and 2CI2 (chymotrypsin
inhibitor 2, I19–83) full unique protein chains. Masks were frozen as ordinal
intervals `[9,13)`, `[23,31)` and `[0,4)` before any model sampling on these
structures. All six source/reference projections passed admission. Original
observations and occupancies remain unchanged; no source repair or trim-to-pass.
HETATM omission is explicit protein-only projection, not whole-entry chemistry.

The first preparation refused genuine observed terminal OXT because its helper
only knew internal canonical atom sets. A separately retained correction admits
**only** an existing true-terminal observed OXT outside the generated mask and
keeps its exact source bytes. No oxygen is deleted, invented or marked generated.
RFD3/Boltz canonical native axes omit this OXT; it is explicitly fixed pass-through
in the final structure, not a model coordinate constraint or native prediction.
Protenix and the MODELLER comparator must represent it. Other missing/unexpected
native atoms fail closed. Do not equate this with identical atom-wise native
conditioning across architectures or qualify strict all-observation inpainting.

These entries are outside the development cases, with distinct family labels.
They are not proven absent from model training. Two entries × three masks do not
make six independent families; five seeds do not make five independent proteins.

## Frozen policy

- Seeds **211, 223, 227, 229, 233**, previously not sampled in this study.
- Models: RFD3, stock Boltz-2 proposal mode, Protenix v1 reinjection, plus the
  unchanged private seeded one-fast-loop MODELLER comparator.
- Full declared protein-chain context and fixed known sequence. No ablation,
  source-coordinate uncertainty weighting, adaptive backend or best-of-five
  candidate selection. Report each backend/seed independently.
- Same pinned checkpoints, owned-source patches and 200 configured timesteps per
  ML backend; RFD3 measured calls are recorded separately. No alternate smaller
  checkpoints, generic MPS CPU fallback or remote MSA/sequence upload.
- Native drift remains separate from restored-source quality. Exact source ATOM
  bytes (including OXT), finite complete requested inventory, unchanged existing
  versioned geometry criteria and global L-chirality determine final acceptance.
- Clean source-restored proposals are not repaired. Failures receive the same
  all-generated-atoms joint CPU1thread repair: max1000 iterations, optimizerseed7,
  no restarts, geometry preconditioning or anchor release. Keep repair errors and
  rejected raw/final structures in denominators; no silent fallback or resampling.
- Held-out coordinates enter post-generation scoring only, never conditioning,
  candidate selection or refinement. Report all-candidate heavy/backbone RMSD
  distributions as well as failures; geometry is not accuracy qualification.
- Positive anchor-motion results cannot change this policy or validation settings.
  Neither successful validation nor a single repaired junction authorizes public
  CLI/Homology/GUI integration, a production-default change or P3.

## Exact budget and execution controls

Six masks × five seeds × three ML models = **90 serial ML calls**; plus
**30 MODELLER** proposals. Existing assets/environments are reused. Each ML call
uses pressure-v2, max1800s, RSS8GiB/MPS6GiB and sustained host-pressure/headroom
guards; CPU jobs use the existing guarded serial executor and study lock.
The worst-case declared per-job wall ceilings sum to 60h; expected execution for
these small shapes is approximately 1–3h, not a guarantee. Estimated new evidence
is below 2GiB; check a 10GiB study-output cap and 2GiB free headroom between jobs.
No new asset download or environment mutation is part of this phase.

Before model weights load, verify the pinned RFD3 builder's exact identity masks,
fixed sequence and zero unconditioned coordinates for all six shapes. Native
transport and every final observation are independently audited. Reserved-seed
entrypoints require the exact `policy.json` SHA256 and explicit phase permission;
old development seed gates remain default-deny for these seeds. All technical
failures stop remaining launches, retaining the unchanged native outputs for an
explicit separately recorded transport recovery if necessary.

`scripts/run_mac_fixed_validation.py` declares immutable policy, jobs, inputs and
source-code bindings before execution. `scripts/mac_validation_inputs.py` owns
reference-blind native identity/frame transport and explicit OXT pass-through.
The numerical inference methods are unchanged; Protenix's private function gains
an explicit allowed-seed parameter with its original development default. No
public API or frozen public request schema is extended.

## Execution lineage

Completed execution evidence root:
`tmp/mac-reconstruction-comparison-2026-11-11/fixed-source-validation-attempt03/`.
Its frozen policy SHA256 is
`8cad45a30ce7e8dc328580bc107a096f140b15751326ff0b544d689c64929434`.
All six inputs passed weightless preflight in all three ML runtimes before this
run launched. Reserved seeds are authorized and execution has started; do not
describe them as untouched. Per-job receipts establish actual seed consumption.

Attempts01/02 are retained **pre-execution** failures, with no checkpoint sampling.
Attempt01 strengthened source/asset binding. Attempt02 exposed two isolated-runtime
transport issues: resolving the venv Python symlink selected the base interpreter,
and upstream Protenix's regular `scripts` package shadowed project namespace imports.
The explicit amendment preserves venv executable paths and imports the policy helper
directly from the declared project `scripts/` path. Numerical methods, scientific
inputs and acceptance thresholds were not changed. Original logs and pre-amendment
entrypoints remain in attempt02. Attempt03 is the first sampling launch, not a
resampling of unsuccessful validation outcomes.

The frozen manifest is an execution declaration, not a success claim. Directory
existence is not completion: require native/execution receipts, per-candidate results,
`ML-completed.json`, the guarded MODELLER completion record and `completed.json`.
An independent auditor must re-read source records (including OXT), inventory,
global chirality, geometry and reference errors before final reporting. Neither
completed calls nor accepted geometry establishes accuracy or public qualification.

Execution and independent readback are now complete: **90/90 ML and 23/30 MODELLER**
final passes. All 120 preserve exact observations/OXT and global L-chirality.
The seven MODELLER rejections remain, including one numerical repair failure; no
resampling or outcome-dependent amendment occurred. Full per-seed results and
hashes are retained at
`tmp/mac-reconstruction-comparison-2026-11-11/fixed-source-validation-independent-readback/`.
See the [comparison log](../research/mac-reconstruction-comparison-2026-11-11.md)
for all-outcome accuracy tables and failure breakdown. This closes the bounded
execution/readback gate, not broad biological or public capability qualification.
