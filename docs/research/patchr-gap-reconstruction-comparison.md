# PATCHR versus dvbfixer gap reconstruction

Date: 2026-10-23. Source-code review, not a reproduced PATCHR benchmark.
PATCHR revision: [`6c161b155b52e1cdee3496b0beac40306793f4c7`](https://github.com/DeepFoldProtein/patchr/tree/6c161b155b52e1cdee3496b0beac40306793f4c7).
The dvbfixer comparison concerns diffusion reconstruction, particularly corrected
Protpardelle-1c cc89; it does not replace the separate Modeller modeling path.

## Main finding

PATCHR's most relevant additional mechanism is **Local Refinement Diffusion
(LRD)**: a second neural denoising pass at template/generated junctions.
dvbfixer currently uses localized **OpenMM refinement**, not neural LRD.
PATCHR's default LRD permits observed boundary atoms to move, unlike dvbfixer's
strict observed-heavy-coordinate preservation. It is therefore not an equivalent
drop-in implementation of the same constrained task.

## How PATCHR constructs missing regions

The template processor builds a full target sequence and maps observed atoms to
it. When applicable, mmCIF `label_seq_id` provides a direct sequence-position
mapping. Otherwise the protein fallback uses local Biopython BLASTP-scored
alignment, with heuristic coverage/match bonuses and author-numbering penalties.
These combined scores are not equivalent to dvbfixer's sequence-first,
peptide-continuity-aware tie breaking. Insertion codes are represented during
input mapping; downstream neural atoms use normalized sequence indices.

A per-atom template mask distinguishes observed from missing atoms, allowing
whole-residue reconstruction and partial-sidechain completion. Boltz-2 and
Protenix are the implemented neural backends; the repository also handles
protein/nucleic-acid/ligand complexes.

Sources: PATCHR
[`alignment.py:409–562`](https://github.com/DeepFoldProtein/patchr/blob/6c161b155b52e1cdee3496b0beac40306793f4c7/scripts/inpainting/processor/alignment.py#L409-L562),
[`template_featurizer.py`](https://github.com/DeepFoldProtein/patchr/blob/6c161b155b52e1cdee3496b0beac40306793f4c7/src/protenix/data/inpainting/template_featurizer.py).
Compare dvbfixer's [shared alignment service](../../src/dvbfixer/sequence_alignment.py).

## Template-Constrained Diffusion (TCD)

PATCHR modifies inference without retraining:

1. Noise template coordinates to the current sampling noise level.
2. Run the neural denoiser.
3. Rigidly align noisy and denoised states using Kabsch.
4. Transform the template with the same rigid transform.
5. After the sampling update, hard-reset masked atoms to the template in the
   evolving coordinate frame.

This is direct coordinate control, not only a template embedding. The noisy
network input and clean post-update state must not be conflated. Nor does
preservation in an evolving rigid frame imply byte-identical final coordinates.

Sources: PATCHR
[`Boltz diffusion`](https://github.com/DeepFoldProtein/patchr/blob/6c161b155b52e1cdee3496b0beac40306793f4c7/src/boltz/model/modules/diffusionv2.py#L572-L954),
[`Protenix sampling`](https://github.com/DeepFoldProtein/patchr/blob/6c161b155b52e1cdee3496b0beac40306793f4c7/src/protenix/model/generator.py#L273-L346).

## Local Refinement Diffusion (LRD)

The default local stage selects a ±3-residue window at template/generated
transitions, starts at sigma 1.5 Å and performs 25 additional denoising steps.
The window includes both observed and generated residues. In the reviewed
implementations, this stage filters to protein atoms and excludes PTMs; general
DNA/RNA support does not imply this same junction refinement is applied to them.
The Boltz implementation also has an internal generated-side-only option,
`boundary_refinement_inpainting_region_mode`, disabled by default.

The Boltz and Protenix local loops are not identical; their shared high-level
name is not proof of numerical equivalence.

Sources: PATCHR
[`Boltz LRD`](https://github.com/DeepFoldProtein/patchr/blob/6c161b155b52e1cdee3496b0beac40306793f4c7/src/boltz/model/modules/diffusionv2.py#L988-L1053),
[`Protenix LRD`](https://github.com/DeepFoldProtein/patchr/blob/6c161b155b52e1cdee3496b0beac40306793f4c7/src/protenix/model/generator.py#L365-L654).

## Observed-coordinate caveat in the reviewed revision

PATCHR's README claims observed coordinates remain exactly as deposited.
The reviewed default code does not establish that guarantee:

- LRD releases template-side boundary atoms explicitly.
- The Boltz final input-frame alignment and exact template overwrite block is
  commented out at lines 1360–1379.
- Protenix defines `_final_align_to_pdb_frame`, but the reviewed sampler return
  paths do not call it.
- A rigid alignment alone cannot undo local observed-atom deformation from LRD.

These are findings about the pinned revision, not claims about every PATCHR
release or a measured displacement in an executed run.

Source: PATCHR
[`commented final restoration`](https://github.com/DeepFoldProtein/patchr/blob/6c161b155b52e1cdee3496b0beac40306793f4c7/src/boltz/model/modules/diffusionv2.py#L1342-L1379).

## Comparison with current dvbfixer

| Concern | PATCHR | dvbfixer |
|---|---|---|
| Main backends in this comparison | Boltz-2 / Protenix | Corrected cc89; also a Protenix adapter |
| Multiple missing regions | Joint represented-complex sampling and boundary refinement | Corrected cc89 samples each gap independently from immutable source-only context, then merges generated coordinates |
| Context | Represented complex | Selected local protein context, at most 512 residues per cc89 engine invocation |
| Per-step coordinate reset | Built into TCD | Profile-specific; implemented for Protenix and as a cc89 experimental callback |
| Corrected public MPS cc89 | Not applicable | Native motif conditioning, `per_step_reinjection=False`, followed by frame synchronization and exact source restoration |
| Junction repair | Neural LRD | OpenMM boundary refinement v8 |
| Observed boundary atoms | Mobile by default during LRD | Kept exact |
| Acceptance | Confidence/diagnostics plus evaluation tools | Independent blocking validation before publication |
| Output metadata | Rebuilt structure, primarily mmCIF, selected template-category merge | Record-preserving PDB materialization with target-bound metadata reconciliation |

dvbfixer already contains [Kabsch plus per-step fixed-coordinate restoration](../../deploy/protenix-v1/reinjection.py).
Do not describe every profile as using it: the
[corrected MPS launcher](../../deploy/protpardelle-1c/mps_runner.py) explicitly
disables cc89 per-step reinjection. Its
[independent-gap orchestration](../../src/dvbfixer/model/diffusion/independent_gaps.py)
does not condition later gaps on earlier generated predictions.

The [OpenMM refiner](../../src/dvbfixer/model/diffusion/boundary_refinement.py)
minimizes generated atoms while freezing observed heavy atoms. The
[independent validator](../../src/dvbfixer/model/diffusion/validate.py) checks
identity/completeness, fixed-coordinate drift, junctions, local lengths/angles,
amide planarity, Ramachandran, sidechain chi1/chi2, severe overlaps, D-Cα,
explicit links and PDB metadata. A runner-success status is not scientific
acceptance.

## Benchmark claims are not a common acceptance rate

PATCHR reports a 940-structure benchmark, missing-region Cα RMSD of 1.781 Å and
99.4% connectivity pass rate. These are author-reported figures from its
[README](https://github.com/DeepFoldProtein/patchr/blob/6c161b155b52e1cdee3496b0beac40306793f4c7/README.md),
not independently reproduced here. The reviewed prediction path does not expose
an equivalent single fail-closed gate matching dvbfixer's full acceptance
contract. Different datasets, allowed observed-atom motion and geometry checks
prevent direct comparison of those percentages with dvbfixer pass rates.

## Implication for the supplied 9CQD failure

The archived refined 9CQD candidate fails at W28→W29, omega −132.2902°:
47.7098° from trans, beyond the unchanged 30° tolerance. Neural local
re-denoising might escape a conformation that local energy minimization does
not fix. It is not established that it will do so with every original atom
frozen. PATCHR's default can also move the observed W29 side of that junction,
which would violate dvbfixer's present coordinate-preservation contract.

See the [9CQD archive audit](multi-gap-diffusion-5j13.md#related-9cqd-metadata-regression-2026-10-23).

## Recommended experiments, not shipped behavior

The requested two-mobility/two-initial-sampling refinement comparison is specified
in the [cc89 LRD experiment plan](../plans/cc89-local-refinement-ablation.md).
The [completed 9CQD results](cc89-lrd-ablation-9cqd.md) have zero accepted
candidates across ten variants, including source-mobile geometry-only checks.
The subsequent generated+immediate-anchor minimization extension passes geometry
in all six additional variants, but all retain strict source-motion failures;
see the same report for the separate policy and measured anchor displacement.
Joint reconstruction of nearby gaps is explicitly deferred there until a
suitable test structure is available.

1. First isolate per-step cc89 reinjection, changing no sequence, context,
   sampling schedule, seed, refinement policy or validation threshold.
2. Then compare generated-only neural LRD, OpenMM-only, and LRD→OpenMM on the
   same raw candidate set. Do not mix initial-sampling changes into this arm.
3. Evaluate all arms with the same independent validator, including exact
   observed-coordinate preservation.
4. Study joint sampling of spatially coupled gaps separately from local repair.

Per-step reinjection is **not** neural LRD, and cc89 reinjection is not a full
reimplementation of PATCHR TCD. Moving observed boundary atoms would require a
separate explicit user policy, not a silent acceptance-threshold change.

### 9CQD seed-7 reinjection ablation

An explicitly private experiment is prepared under the ignored directory
`test/errors_diffusion/9cqd-reinjection-seed7/`. It retains the archived request
and source byte-for-byte, verified against the archive checksum manifest:
targets M/W/X 27/214/221; gaps W26–28 and X140–141; seed 7. It retains
500 steps, step scale 1.2, s_churn 200, neutral chain phase, physical position
features, independent source-only gap contexts and isolated CPU refinement v8.

The production profile remains unchanged. A saved private adapter copy removes
the reinjection prohibition and allows the explicitly declared source-only
context policy in the research request; `adapter.diff` and its digest record
these changes. The execution request changes only the profile label to the
existing source-physical research identity: sequences, masks, gaps and context
plans are asserted equal to the archived request. This is necessary because
the public-profile trace contract correctly forbids reinjection. The contract
itself is not changed.
The frozen Apple engine does not include the callback hook, so a private engine
copy additionally receives the repository's existing `per-step-callback.patch`.
The private environment verifier checks the combined Apple+callback model-file
digest instead of the production Apple-only digest. The manifest records both
the callback patch and combined engine hashes. Existing callback code must
demonstrate exactly 500 updates per gap and zero post-projection error in the
sampler's motif frame. Provenance labels the run as a private reinjection
ablation, not a qualified public profile.
Current metadata reconciliation applies to the new output and the baseline
metadata audit. The original archive is never overwritten.

The common aggregator also intentionally rejects research-labelled multi-gap
traces. Its already-written merged coordinates and independently checked
per-gap traces are retained. Scientific refinement is called separately with
the same CPU implementation and parameters as the production wrapper, without
forging an aggregate production trace or publishing a production bundle.
Execution and scientific acceptance results are recorded separately when
available; setup alone is not evidence of a successful reconstruction.

Pre-sampling attempts retained in the experiment directory stopped first on
cached bytecode detected by the strict source-tree verifier, then on the missing
`step_callback` sampler argument. Neither attempt produced a reconstruction;
these checks were not bypassed in the production runtime.
An additional attempt completed gap-0 sampling but was rejected while building
a trace with the incompatible public profile label. That evidence is retained;
the final experiment uses the research label rather than weakening the contract.

### Completed result: reinjection did not resolve this candidate

Both gaps completed on MPS with **500 callbacks each**, zero maximum
post-projection error in the motif frame, and 1,000 observed first-layer RoPE
calls per gap. Sampling took approximately 19.3 and 18.3 seconds. CPU refinement
v8 took approximately 30.0 seconds, with the same 500-iteration limit, eight
restarts, 0.25 Å perturbation and seed 7 as the archived short-gap run.

Independent scientific validation, after normalizing only archived metadata for
a fair geometry comparison:

| Metric | Archived refined, no reinjection | Reinjection, refined |
|---|---:|---:|
| Scientific validation passed | No | No |
| Failing gate | Amide planarity | Amide planarity |
| Non-planar amides near generated regions | 1 | 2 |
| W28→W29 omega | −132.2902° | −116.0406° |
| W28→W29 deviation from trans | 47.7098° | 63.9594° |
| Junction C–N minimum/maximum | 1.2258 / 1.3433 Å | 1.2317 / 1.3451 Å |
| Severe overlaps | 0 | 0 |
| Detectable D-Cα | 0 | 0 |
| Maximum observed-atom displacement | 0 Å | 0 Å |

The second non-planar amide in the reinjection candidate is W25→W26,
omega 146.2973° (33.7027° from trans). The unchanged tolerance is 30°.
Both raw arms fail the same three gates: junction peptide connectivity, amide
planarity and severe overlap. Raw W28→W29 omega is −123.3654° in the archive
versus −123.0896° with reinjection. Thus final geometry was not rescued simply
by projecting observed atoms at every step. One archived-control comparison
does not establish population-level superiority/inferiority; MPS sampling and
CPU minimization are not assumed bit-identical across executions or machines.

The rejected research output is retained at
`test/errors_diffusion/9cqd-reinjection-seed7/refined/candidate.pdb`, with
`comparison.json`, `refinement-report.json`, per-gap traces, source/adapter
copies and failed-attempt logs beside it. It was **not published as an accepted
model**. No validation threshold or production profile was changed.

Evidence digests:

| Artifact | SHA-256 |
|---|---|
| Original archived request | `713d5fbc971b46bdc69cbf8a4215e216d8747b1cedf39b7d6f8e464c33e32246` |
| Research-labelled execution request | `f0d3f96557f3223a432d1c18a722792393fce5ca1d0f5604b3aa3d1f5b66b9d4` |
| Merged reinjection raw PDB | `7453570841f46986a782d74501c730cdb01c56f34a5ebb190aa0e00745fd95a5` |
| Reinjection refined PDB | `94f4813d3c7b5eb26fc9ea44122ae50fa6acada9195f2ab142869a7ebb46fca7` |
| Private Apple+callback models.py | `5dc3ba8dd6a0bd77c9c76639fad05cf227d5ff4ae2bd490a1f272430bccaa66a` |

The production Apple models.py remains at its frozen
`8513fac0d18d3080b14673bb9b6ed3b231e1d31d30449d3b8c6cc09f0445800b`
digest. Callback/adapter tests: 52 passed. Documentation local links and agent
documentation maps validate successfully.
