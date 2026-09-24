# Diffusion Gap-Reconstruction Implementation Plan

- Status: proposed.
- Branch: `diffusion`.
- Scope owner: Structure Preparation bounded context.
- Change group: `missing-atom-rebuild-chirality`.
- Coordination rule: assign one exclusive coordinating owner before implementation.
- Related research:
  - [`../research/reconstruction-and-modeling-backends.md`](../research/reconstruction-and-modeling-backends.md)
- Production baselines retained:
  - Salilab MODELLER for `model` and `homology`.
  - PDBFixer/OpenMM Modeller in the documented legacy preparation path.
  - `tleap-reduce` for its documented pure-protein preparation scope.

## Goals

- [x] Add an internal experimental, backend-neutral diffusion protocol path for sequence-guided protein gap reconstruction; no public engine backend is claimed.
- [x] Preserve deposited coordinates outside explicitly generated or movable regions through exact identity and fixed-heavy coordinate gates.
- [ ] Implement the PATCHR-like method independently rather than importing, wrapping, vendoring, or copying PATCHR. CPU synchronization/reinjection primitives and conformance gates exist; real per-step sampler integration remains externally blocked.
- [x] Keep ML frameworks, CUDA libraries, and model-specific packages out of the core DVBFixer environment.
- [x] Record enough available provenance to reproduce or audit every generated candidate without inventing missing engine evidence.
- [x] Validate identity, connectivity, supported geometry, chirality, and clashes before publishing a candidate.
- [ ] Evaluate mosaic-first diffusion for comparative/homology modeling only after local gap reconstruction passes its acceptance gates.

## Non-goals

- [x] Do not replace or deprecate MODELLER in this work.
- [x] Do not change the default `model` or `homology` backend during the research phases.
- [x] Do not add automatic fallback from diffusion to MODELLER.
- [x] Do not use full-chain prediction as a substitute for fixed-coordinate gap repair.
- [x] Do not treat missing heavy atoms inside an otherwise present residue as the same problem as missing whole residues.
- [x] Do not move CIF parsing into the scientific modeling stage.
- [x] Do not claim support for ligands, glycans, PTMs, cofactors, metals, or covalent non-peptide links until each class has separate evidence.
- [x] Do not store model weights in Git.

## Initial Supported Slice

- [x] Accept canonical L-protein inputs only at the internal scope boundary.
- [x] Require a complete, known target sequence in the request; FASTA/SEQRES extraction remains owned by the existing CLI/model boundary.
- [x] Require unambiguous shared sequence placement.
- [x] Support internal gaps with two observed anchors.
- [x] Start with gap lengths of 3-12 residues.
- [x] Reject terminal one-anchor gaps in the first slice.
- [x] Reject alternate-location ambiguity in or adjacent to the generated region.
- [x] Reject multiple `MODEL` blocks in the first slice.
- [x] Reject retained unsupported heterogens or external covalent links near the generated region.
- [x] Return an explicit `unsupported` result without creating a public PDB for every out-of-scope case. The internal pipeline returns a typed unsupported outcome before runner launch and publication tests prove no destination bundle or staging orphan is created.

## Mandatory Identity And Publication Rules

- [x] Preserve case-sensitive chain IDs.
- [x] Preserve residue insertion codes.
- [x] Use `(chain, resid, icode, atom)` as atom identity.
- [x] Use `(chain, resid, icode)` as residue identity.
- [x] Preserve all atoms and explicit links outside the generated mask.
- [x] Preserve `keepIds=True` on every OpenMM `PDBFile.writeFile` call. The diffusion core does not currently serialize through `PDBFile.writeFile`; the repository-wide rule remains unchanged.
- [x] Normalize CIF input only at the CLI boundary. The diffusion core accepts only normalized PDB bytes and adds no CIF reader.
- [x] Use `dvbfixer.ffutils.dat.DatRecord` for the downstream `.dat` sidecar.
- [x] Publish a candidate only after every hard validation gate passes.
- [x] Write candidate PDB, `.dat`, diffusion provenance, and a bundle index atomically as one logical result through one directory rename.
- [x] Leave no partial public output after runner, parsing, validation, or publication failure.

## Engine And License Decisions

### RFdiffusion v1 Baseline

- [x] Select RFdiffusion v1 as the first external backbone-inpainting benchmark adapter. The internal adapter is implemented and GPU-smoked; it remains outside public dispatch.
- [x] Pin source revision `bf42b54c20a99dd7350456c85985ed4d83b95d48` in the research inventory.
- [x] Pin the exact checkpoint to independently measured SHA-256 `0fcf7d7c32b4848030aca3a051e6768de194616f96ba6c38186351a33bfc6eca`. Upstream still provides no signed checksum manifest, so publisher-authenticated supply-chain verification remains unavailable.
- [x] Archive BSD-3-Clause evidence for code and upstream license commit `820bfdfaded8c260b962dc40a3171eae316b6ce0`, which explicitly covers README-linked weights. A final redistribution review remains required before embedding weights in a service image.
- [x] Treat any future RFdiffusion output as a backbone-level benchmark until side-chain materialization and all-heavy-atom validation are complete.
- [x] Do not expose raw RFdiffusion output as a production DVBFixer repair result.

### All-Atom Feasibility Spike

- [ ] Evaluate Protenix v1 first; revision `85767b811c40ed46e73a9b39519cf6bfca8701ba` and its Apache-2.0 source claim are inventoried, but checkpoint acquisition and sampler-hook evidence remain blocked.
  - [x] Pin the audited v1 source revision and record the upstream Apache-2.0 code/model-parameter claim; redistribution review is still required.
  - [x] Do not substitute Protenix v2 without a new license review.
  - [ ] Verify whether the denoising state can be intercepted before every next step.
- [ ] Evaluate Boltz-2 if Protenix v1 does not provide a maintainable hook; revision `b1ebfc46ecf57f5414e0d1a6f9027bbb122c53bc` is inventoried, but checkpoint and sampler-hook evidence remain blocked.
  - [x] Pin the audited MIT source revision and record the upstream code/weight claim; exact weight revision, URL, hash, and redistribution review remain unresolved.
  - [x] Treat ordinary template conditioning as insufficient for exact fixed-coordinate preservation.
  - [ ] Verify whether fixed-coordinate reinjection can be implemented without an unmaintainable upstream fork.
- [x] Add a backend-neutral conformance record that refuses to claim reinjection unless mutable state, stable identity mapping, and exact fixed-coordinate overwrite are exposed at every denoising step.
- [x] Stop before user-facing integration if neither pinned engine provides stable denoising hooks.
- [ ] Record whether a maintained fork, upstream API change, or different sampler would be required after the hook spikes run.

### PATCHR Comparator

- [x] Use PATCHR only as an external scientific comparator.
- [x] Do not import PATCHR into DVBFixer.
- [x] Do not vendor PATCHR source.
- [x] Do not wrap the PATCHR CLI as the DVBFixer implementation.
- [x] Require a separate ADR before changing this repository policy.

### Excluded Distribution Candidates

- [x] Exclude FrameDiPT from the distributable backend because of CC BY-NC-SA terms.
- [x] Exclude Chroma parameters from the distributable backend because of non-commercial parameter terms and API-key distribution.
- [x] Exclude RFdiffusion2 until official checkpoint rights and optional Chai dependencies are unambiguous.
- [x] Exclude engines with unclear checkpoint redistribution rights from the shipped adapter list.
- [x] Keep unconditional-generation or motif-scaffolding engines as research comparators unless they demonstrate the required two-anchor repair contract.

### Per-Artifact License Inventory

- [x] Record code repository and revision.
- [x] Record source-code license and exact license-text revision.
- [x] Record checkpoint URL and independently measured SHA-256; retain the absence of a publisher-signed checksum as explicit provenance.
- [ ] Record checkpoint license and redistribution decision. Upstream claims and unresolved review status are recorded; no engine is approved for distribution.
- [ ] Record auxiliary model and cache URLs, hashes, and licenses. Required auxiliary artifacts remain unresolved.
- [ ] Record container base-image and package-lock provenance. Required digests and lock hashes remain unresolved.
- [x] Record relevant training-data cutoff claims for leakage analysis, including unresolved sequence-level membership.

## Backend-Neutral Contract

- [x] Add `src/dvbfixer/model/diffusion/contract.py`.
- [x] Define typed `AtomIdentity` and `ResidueIdentity` values.
- [x] Define a `GapRegion` with:
  - [x] target-chain identity;
  - [x] target-sequence interval;
  - [x] left and right anchors;
  - [x] generated residues;
  - [x] movable junction window.
- [x] Define a versioned `DiffusionRequest` with:
  - [x] normalized PDB input;
  - [x] target sequences;
  - [x] shared sequence placement;
  - [x] explicit gap ranges;
  - [x] fixed and generated atom masks;
  - [x] retained explicit links;
  - [x] candidate count;
  - [x] seed list;
  - [x] backend options.
- [x] Define a `DiffusionCandidate` with:
  - [x] candidate coordinate artifact;
  - [x] generated atom and residue identities;
  - [x] raw backend score;
  - [x] score provenance;
  - [x] warnings.
- [x] Define a `DiffusionResult` with:
  - [x] `success`, `unsupported`, or `failed` status;
  - [x] candidates;
  - [x] validation summaries;
  - [x] runner diagnostics;
  - [x] backend provenance.
- [x] Add strict JSON serialization and schema-version validation.
- [x] Reject unknown or incompatible schema versions explicitly.

## Mask And Sequence Placement Layer

- [x] Add `src/dvbfixer/model/diffusion/masks.py`.
- [x] Reuse the existing shared sequence-placement behavior instead of implementing a second alignment policy.
- [x] Construct fixed/generated masks from the authoritative alignment.
- [x] Preserve insertion-code-aware residue identity.
- [x] Preserve case-distinct chains such as `D` and `d`.
- [x] Define deterministic left and right anchor windows.
- [x] Define a separately bounded movable junction window.
- [x] Reject ambiguous placements instead of selecting a scientifically silent alternative.
- [x] Reject one-anchor and no-anchor regions in the initial slice.

## External Runner Boundary

- [x] Add `src/dvbfixer/model/diffusion/runner.py`.
- [x] Use a versioned subprocess protocol:
  - [x] DVBFixer writes `request.json` into an isolated workspace.
  - [x] The external runner writes a raw `RunnerResult` in `result.json` plus candidate artifacts; independent `DiffusionResult` validation summaries remain DVBFixer-owned.
  - [x] DVBFixer validates every returned path and digest.
- [x] Enforce workspace path containment.
- [x] Reject symlink and path-traversal escapes.
- [x] Bound runtime and captured output.
- [x] Capture exit code, stdout, and stderr without leaking credentials.
- [x] Treat malformed or missing result manifests as failures.
- [x] Add typed, CPU-testable local preflight results for missing runner, immutable image identity, local checkpoint, checkpoint digest mismatch, CUDA device, and incompatible protocol version. Engine-specific launchers still must supply real image/device discovery in Phase 2.
- [x] Add a deterministic fake CPU runner for unit and integration tests.
- [x] Keep runner-specific Python packages outside the core DVBFixer dependency set.

## PATCHR-Like Constrained Sampling

- [ ] Build the complete target topology before sampling.
- [x] Separate fixed observed atoms from generated missing-region atoms in the backend-neutral request/mask contract.
- [ ] Initialize or sample only the generated region.
- [ ] After every denoising update:
  - [x] provide a CPU-tested weighted Kabsch synchronization primitive from reliable fixed anchors;
  - [x] provide a CPU-tested generated-frame synchronization primitive;
  - [x] provide a CPU-tested exact fixed-coordinate overwrite primitive;
  - [x] verify atom identity equality in the synchronization/reinjection primitive;
  - [ ] integrate all four operations inside every update of a pinned real sampler.
- [x] Run a localized post-sampling OpenMM pass around both peptide junctions for the RFdiffusion baseline; a real-sampler second diffusion pass remains pending.
- [x] Keep atoms outside the generated residue set fixed during the localized pass. Fixed request atoms take precedence over the broader movable-junction support window.
- [x] Materialize all expected canonical heavy atoms in generated residues for the RFdiffusion baseline.
- [x] Leave hydrogen placement and protonation to the downstream preparation stage contract.
- [x] Define three explicit ablation modes and their required capability evidence:
  - [x] template conditioning only;
  - [x] template conditioning plus per-step reinjection;
  - [x] reinjection plus local boundary refinement.
- [ ] Execute those ablations with a pinned real sampler.

## Independent Validation And Ranking

- [x] Add `src/dvbfixer/model/diffusion/validate.py`.
- [x] Compare atom and residue identities outside the generated mask exactly.
- [x] Compare explicit links outside and across the generated region.
- [x] Measure fixed-heavy-atom RMSD and maximum displacement.
- [x] Validate peptide connectivity at both junctions.
- [ ] Validate generated bond lengths, angles, planarity, and atom completeness. Canonical heavy-atom completeness, generated C-N break gates, severe local bond-length outliers, broad canonical backbone/peptide-angle gates, and peptide-amide planarity are implemented through shared diagnose policies; finer residue-specific side-chain geometry remains pending.
- [x] Validate severe intra-region and region-context clashes.
- [ ] Validate Ramachandran and rotamer quality. Gross general-residue Ramachandran outliers are rejected through the bundled MDAnalysis Lovell-reference 99% contour, and gross pooled χ1/χ2 outliers are rejected through the bundled Janin 98% contour with conservative periodic neighborhoods; GLY/PRO/pre-PRO Ramachandran, χ1-only residues, χ3-χ5, and residue/backbone-dependent full rotamer quality remain pending.
- [ ] Run `fix_ca_chirality` only as an explicitly recorded repair step. Validation currently rejects D geometry and performs no silent repair.
- [x] Run `assert_all_l` after the final heavy-atom coordinate change.
- [x] Reject any candidate that fails a hard gate.
- [x] Rank only passing candidates.
- [x] Keep backend confidence, closure metrics, geometry metrics, and physical energies as separately named score components. First-slice ranking keeps DVBFixer geometry metrics separate and does not use raw backend confidence as a shared-scale metric.
- [x] Do not compare backend confidence directly with MODELLER `molpdf` as though they share a scale.

## Provenance Manifest

- [x] Write `<stem>.diffusion.json` beside each published diffusion candidate inside its committed result bundle.
- [x] Include a manifest schema version.
- [x] Record DVBFixer version and commit, using explicit `unknown` outside a verifiable Git checkout.
- [x] Record runner protocol version.
- [x] Record engine repository and revision.
- [x] Record optional container digest or environment-lock identity/hash without inventing absent values.
- [x] Record optional checkpoint hash and source/checkpoint license identifiers without inventing absent values; auxiliary-artifact inventory remains adapter-specific.
- [x] Record request hash, target sequence, gap masks, and fixed/generated identity sets.
- [x] Record seed and optional device, precision, framework, CUDA, driver, and deterministic flags.
- [x] Record known nondeterministic kernels or operations when supplied by the runner.
- [x] Record raw score and score provenance.
- [x] Record every validation metric and pass/fail decision.
- [x] Record final PDB and `.dat` SHA-256 values.
- [x] Keep diffusion environment metadata out of `DatRecord` unless a later general sidecar decision explicitly changes that contract.

## Atomic Publication

- [x] Add `src/dvbfixer/model/diffusion/pipeline.py`.
- [x] Generate and validate candidates entirely inside a temporary workspace.
- [x] Select the ranked passing candidates before touching public destinations.
- [x] Build candidate-matched `.dat` records through `DatRecord`.
- [x] Stage PDB, `.dat`, provenance manifest, and bundle index in a hidden same-parent directory.
- [x] Commit the complete candidate artifact set atomically with one no-replace directory rename; a destination created concurrently is preserved and publication fails closed.
- [x] Remove staged or post-rename artifacts after any failed commit/durability step.
- [x] Preserve the source input on every failure path.

## Implementation Phases

### Phase 0: Policy, Corpus, And Thresholds

- [x] Update the research note with engine, licensing, and hardware findings.
- [x] Add or update the DDD implementation task without presenting the backend as shipped.
- [x] Accept an ADR covering:
  - [x] subprocess/container isolation;
  - [x] separate diffusion provenance;
  - [x] no automatic scientific fallback;
  - [x] MODELLER remaining the default;
  - [x] independent PATCHR-like implementation.
- [x] Select reviewed benchmark fixtures.
- [x] Record fixture provenance and checksums.
- [x] Record exact masks and expected supported/unsupported classification.
- [x] Record upstream training cutoffs and sequence-identity leakage metadata where available.
- [x] Freeze numerical thresholds before comparative inference runs.

### Phase 1: CPU-Testable Core

- [x] Implement the versioned request/result contract.
- [x] Implement deterministic mask construction.
- [x] Implement weighted Kabsch and coordinate-reinjection primitives.
- [x] Implement the isolated runner protocol.
- [x] Implement the deterministic fake runner.
- [x] Implement independent validation and ranking for the accepted initial hard gates. Finer residue-specific geometry, GLY/PRO/pre-PRO Ramachandran, χ1-only and χ3-χ5 torsions, and residue/backbone-dependent full rotamer validation remain explicitly outside the current claim.
- [x] Implement atomic publication.
- [x] Implement separate provenance manifests.
- [x] Keep public CLI behavior unchanged through the runner slice.
- [x] Keep MODELLER execution unchanged through the runner slice.

### Phase 2: RFdiffusion v1 GPU Baseline

- [x] Build a pinned native Linux/NVIDIA runner environment for the first benchmark host; immutable service-container packaging remains separate.
- [x] Add checksum-verified checkpoint acquisition instructions.
- [x] Map the supported one-chain/one-gap `DiffusionRequest` slice to RFdiffusion contig/inpainting inputs.
- [x] Map generated coordinates back to stable DVBFixer identities through target-sequence ordinals.
- [x] Materialize canonical side chains through the explicitly recorded PDBFixer pure-protein route.
- [x] Re-run all all-heavy-atom validation after side-chain materialization.
- [x] Keep outputs inside benchmark workspaces until hard gates pass.
- [x] Measure fixed-atom drift, closure, withheld quality, runtime, RAM, VRAM, same-seed raw/refined repeatability, and seed variability for the reviewed 5- and 10-residue smoke cases. Raw candidates failed junction/chirality/clash gates; repeated v3 post-refinement candidates pass all hard gates and are byte-identical within each same-seed pair.
- [ ] Build and verify an immutable NVIDIA Docker runner image after selecting a base-image digest. Keep the checkpoint mounted and digest-verified by default until redistribution review explicitly permits embedding it.

### Phase 3: All-Atom Constrained Sampler

- [ ] Complete the Protenix v1 hook feasibility spike.
- [ ] Complete the Boltz-2 hook spike only if needed.
- [ ] Select a sampler only after per-step control is demonstrated.
- [x] Implement the backend-neutral weighted Kabsch synchronization primitive; real-sampler per-step integration remains pending.
- [x] Implement the backend-neutral exact fixed-coordinate reinjection primitive; real-sampler per-step integration remains pending.
- [x] Implement localized post-sampling boundary refinement for the RFdiffusion baseline; this does not satisfy the per-step all-atom sampler ablation.
- [ ] Run the three-way ablation benchmark.
- [ ] Require complete canonical heavy atoms.
- [ ] Require final `assert_all_l` after every candidate's last heavy-coordinate change.

### Phase 4: Experimental `model` CLI

- [ ] Add `--backend {modeller,diffusion}` only after Phase 3 hard gates pass.
- [x] Keep `modeller` as the default; no public diffusion dispatch currently exists.
- [ ] Add diffusion runner/checkpoint/device options in a separate argparse group.
- [ ] Reject incompatible existing options before launching the runner.
- [ ] Preserve FASTA/SEQRES and content-selection semantics where the diffusion scope supports them.
- [x] Do not fall back automatically to MODELLER.
- [ ] Add backend preflight to `doctor` without removing its stable report sections.
- [ ] Regenerate CLI reference and GUI command schema through their generators.
- [ ] Do not propagate diffusion options through `zbs` in this phase.
- [ ] Expose diffusion through the existing asynchronous managed-job API only after the experimental CLI gates pass; do not add a synchronous inference route or duplicate adapter science in Node.

### Phase 5: Mosaic-First Homology

- [ ] Keep `selected_template_mosaic.pdb` as the authoritative coordinate frame.
- [ ] Export companion coverage metadata from `materialize_template_plan` instead of recomputing template ownership later.
- [ ] Preserve zero-based half-open template-plan masks.
- [ ] Mark covered template atoms as fixed.
- [ ] Generate only uncovered insertions, substitutions, and bounded junction windows.
- [ ] Start with one-chain internal insertions.
- [ ] Add multi-chain and multi-template cases only after one-chain acceptance.
- [ ] Preserve distinct antibody H/L chains and insertion codes.
- [ ] Compare against MODELLER using the same target and template plan.
- [ ] Use OpenFold or Boltz full-chain prediction only as an independent plausibility comparator.
- [ ] Do not expose a public homology diffusion backend until mosaic adherence and multi-chain gates pass.

### Phase 6: Production Decision

- [ ] Complete code, weight, auxiliary-artifact, and container license review.
- [ ] Complete representative GPU deployment and reliability evaluation.
- [ ] Compare against MODELLER across every claimed supported stratum.
- [ ] Document supported and unsupported inputs.
- [ ] Make production promotion a separate evidence-backed decision.
- [x] Do not remove MODELLER or PDBFixer as a side effect of any future promotion.

## Benchmark Corpus

- [x] Declare a reviewed withheld-coordinate internal 3-5-residue case.
- [x] Declare a reviewed withheld-coordinate internal 6-12-residue case.
- [ ] Add a separate research stratum for gaps of 13-25 residues.
- [x] Expand beyond the initial 8CZ8 cases to an independent 8B01 regular-loop stratum with reviewed 5-/10-residue masks, same-seed repeats, independent hard-gate validation, and a same-configuration MODELLER comparator.
- [x] Include isolated target-chain-only glycine/proline-rich 7X35 masks (`AGQGP` and `PPGGPVP`) with same-seed repeats, independent validation, and MODELLER comparators. The partner chains and PLM are explicitly excluded, so these do not satisfy interface/ligand-context coverage.
- [x] Declare a reviewed interface-adjacent 7X35 `AWVPR` mask with retained
  protein partner chain `B`, a `2.36 Å` generated-region partner contact, and no
  retained PLM. The adapter now supplies that partner as fixed RFdiffusion
  receptor context and fail-closes on post-sampling partner-backbone drift.
- [ ] Complete pinned A100 same-seed, independent-validation, and MODELLER
  measurements for the declared interface-adjacent case.
- [ ] Include antibody insertion codes.
- [x] Exercise case-sensitive chain identity in CPU unit tests; add a reviewed benchmark structure with case-distinct chains before real-engine claims.
- [x] Declare terminal one-anchor gaps only as a separate unsupported/later stratum.
- [x] Declare retained-heterogen cases as a separate unsupported/later stratum; add covalent-link fixtures when that stratum is reviewed.
- [x] Keep current reviewed structures under `tests/fixtures/`.
- [x] Document current fixture provenance in `tests/fixtures/README.md`.
- [x] Validate current fixture checksums against `tests/fixtures/MANIFEST.sha256`; regenerate it whenever fixtures change.

## Initial Hard Gates

- [x] Preserve 100% of atom and residue identities outside the generated mask.
- [x] Preserve 100% of applicable explicit links.
- [x] Require fixed-heavy-atom RMSD no greater than `0.01 Å` after final serialization.
- [x] Require fixed-heavy-atom maximum displacement no greater than `0.03 Å`.
- [x] Require zero detectable D-C-alpha centers.
- [x] Require complete expected heavy atoms for every generated canonical residue.
- [x] Require both junction C-N distances to fall within `1.20-1.45 Å`.
- [x] Reject any backbone break greater than `1.8 Å` in the generated region.
- [x] Reject severe steric overlaps in the generated and junction neighborhoods.
- [x] Reject gross general-residue Ramachandran outliers in the generated and junction neighborhoods; class-specific GLY/PRO/pre-PRO validation remains pending.
- [x] Reject gross pooled χ1/χ2 side-chain torsion outliers where both torsions are defined; this is not complete residue/backbone-dependent rotamer validation.
- [x] Require unsupported or ambiguous inputs to create no public PDB or result bundle.
- [x] Provide a CPU-testable same-seed repeatability classifier at `0.01 Å` and apply it to independent pinned A100 workspaces. Both initial v3 refined mask pairs measure `0.0 Å` generated-heavy RMSD and maximum displacement.
- [x] Version threshold changes instead of adjusting them after viewing benchmark outcomes. The current threshold is a named module constant and recorded in each repeatability assessment.

## Comparative Metrics

- [x] Compute gap backbone RMSD against withheld coordinates in the benchmark API and record pinned real-engine values for the initial 5-/10-residue subset.
- [x] Compute gap all-heavy-atom RMSD in the benchmark API and record pinned real-engine values for the initial 5-/10-residue subset.
- [ ] Report lDDT and GDT-HA or TM-score where meaningful after a reviewed implementation and real candidates exist.
- [x] Compute fixed-heavy and anchor-heavy RMSD.
- [x] Summarize peptide closure and independent-validation pass rates from explicit per-candidate decisions.
- [x] Preserve bond, angle, planarity, Ramachandran, pooled χ1/χ2, chirality, and clash metrics from independent validation manifests.
- [x] Compute top-1 and oracle top-k backbone quality from a declared ranking order.
- [x] Compute pairwise candidate diversity and ranking enrichment.
- [x] Compute success, unsupported, and failure rates by stratum from explicit run outcomes; no real-engine rates are claimed yet.
- [x] Aggregate observed wall time, model-load time, peak RAM, and peak VRAM while preserving missing values instead of fabricating them; no real-engine values are claimed yet.
- [x] Compute external-process timeout and crash rates from explicit failed-run evidence; no real-engine rates are claimed yet.
- [ ] Report conditioning/reinjection/boundary-refinement ablation results.
- [x] Encode and measure the gap-backbone RMSD gate as no more than `0.25 Å` worse than MODELLER. Refined RFdiffusion beats the MODELLER median on all six measured 8CZ8/8B01/7X35 masks; MODELLER's best `PPGGPVP` candidate is nevertheless better than the single RF candidate (`3.54 Å` versus `4.28 Å`).
- [x] Encode and measure the junction-pass gate as no lower than MODELLER. Both measured backends pass both junctions on all six masks after refinement.
- [x] Encode and measure fixed-coordinate adherence as strictly better than MODELLER: refined RFdiffusion is `0.0 Å` on all six masks, while every measured MODELLER comparator moves deposited coordinates.

## Hardware And CI Matrix

### CPU-Only Linux

- [x] Run contract serialization tests in the focused diffusion suite.
- [x] Run mask and identity tests in the focused diffusion suite.
- [x] Run Kabsch and reinjection primitive tests in the focused diffusion suite.
- [x] Run fake-runner tests, including same-seed measurement across independent workspaces.
- [x] Run provenance, checksum, local immutable-artifact preflight, and license-inventory tests.
- [x] Run timeout, crash, containment, hard-link/symlink, and atomic no-replace publication tests.
- [x] Run PDB validation, `.dat`, chirality, and benchmark-scoring tests.
- [x] Treat any tiny CPU inference as an optional smoke test only when upstream officially supports it; no real-engine CPU inference is claimed.
- [x] Do not use CPU inference for representative ensembles or acceptance benchmarks.

### macOS Apple Silicon — Deferred Follow-Up

Native Apple Silicon diffusion inference is intentionally deferred until the
Linux/NVIDIA implementation has a stable main framework and has completed its
Phase 2 backbone baseline and Phase 3 all-atom constrained-sampler tests. It is
not part of the current Phase 2-4 implementation cycle and must not delay or
weaken the A100 acceptance gates.

- [ ] Reassess pinned-engine CPU/MPS feasibility only after the main framework and Linux/NVIDIA Phase 2/3 tests are stable.
- [ ] Run the same core, fake-runner, and validation suites as CPU Linux on an actual Apple Silicon host before claiming platform support.
- [ ] Define a separate Apple Silicon research matrix covering upstream arm64 packages, unsupported operators, precision, unified-memory use, repeatability, and runtime.
- [x] Treat MPS inference as exploratory unless the selected pinned engine officially supports it.
- [x] Do not use MPS output as a release acceptance gate or as evidence for the Linux/NVIDIA gates.
- [x] Keep macOS usable as a DVBFixer core/client platform; a future remote Linux/NVIDIA runner may remain the practical inference path.
- [x] Document that Docker Desktop on macOS does not provide NVIDIA CUDA.

### Linux/NVIDIA

- [x] Use Linux x86_64 with NVIDIA driver 535.104.05 and the pinned native CUDA 11.1 environment; the future container runtime is not yet pinned.
- [x] Start short-gap evaluation on an A100-SXM4-40GB.
- [ ] Prefer 48-80 GB VRAM for larger complexes, all-atom models, multi-sample ensembles, and profiling.
- [x] Provide at least 64 GB system RAM and fast local temporary/model-cache storage.
- [x] Record actual peak memory instead of treating 24 GB as an upstream guarantee. The first instrumented run observed approximately 5.96 GB RAM and 3.70 GB VRAM.
- [x] Add a documented representative smoke profile before selecting long-term hardware.

### CI Separation

- [ ] Make the CPU contract lane mandatory in repository CI configuration; the suite itself is CPU-only and passing locally.
- [x] Keep the GPU inference lane opt-in or self-hosted until its environment is stable; no required GPU lane exists.
- [x] Return a clear typed preflight issue when no CUDA device is reported by an adapter launcher.
- [x] Do not let missing CUDA fail the core test suite.

## Expected Files

- [x] Add `src/dvbfixer/model/diffusion/__init__.py`.
- [x] Add `src/dvbfixer/model/diffusion/contract.py`.
- [x] Add `src/dvbfixer/model/diffusion/masks.py`.
- [x] Add `src/dvbfixer/model/diffusion/scope.py`.
- [x] Add `src/dvbfixer/model/diffusion/geometry.py`.
- [x] Add `src/dvbfixer/model/diffusion/runner.py`.
- [x] Add `src/dvbfixer/model/diffusion/validate.py`.
- [x] Add `src/dvbfixer/model/diffusion/pipeline.py`.
- [x] Add `src/dvbfixer/model/diffusion/provenance.py`.
- [x] Add `src/dvbfixer/model/diffusion/benchmark.py`.
- [x] Add `src/dvbfixer/model/diffusion/preflight.py`.
- [x] Add `src/dvbfixer/model/diffusion/sampler.py`.
- [ ] Change `src/dvbfixer/model/cli.py` only in the experimental CLI phase.
- [ ] Change `src/dvbfixer/model/pipeline.py` only for explicit backend dispatch after gates pass.
- [ ] Change `src/dvbfixer/doctor.py` for runner, checkpoint, and device preflight.
- [ ] Change `src/dvbfixer/homology_plan.py` only in the mosaic-first phase.
- [ ] Change `src/dvbfixer/homology.py` only after the homology contract passes focused tests.
- [ ] Avoid changing `src/dvbfixer/ffutils/dat.py` unless a general sidecar requirement is demonstrated.
- [x] Add the isolated native runner environment lock and license inventory; an immutable Docker base digest remains pending before service-image acceptance.
- [x] Add focused diffusion unit and integration tests for the contract, masks, geometry, research inventory, runner protocol, fake runner, independent validation, publication, provenance, adapter preflight, sampler conformance, repeatability, and benchmark metrics.
- [x] Add reproducible internal builders/analyzers for frozen withheld-coordinate requests, same-seed workspace pairs, and MODELLER comparator outputs.
- [ ] Add reviewed fixtures and regenerate their manifest when benchmark structures are added.

## Focused Verification

- [x] Run the diffusion core tests:

  ```bash
  pytest -q tests/test_diffusion_benchmark.py tests/test_diffusion_contract.py \
    tests/test_diffusion_geometry.py tests/test_diffusion_masks.py \
    tests/test_diffusion_pipeline.py tests/test_diffusion_preflight.py \
    tests/test_diffusion_provenance.py \
    tests/test_diffusion_research_inventory.py tests/test_diffusion_runner.py \
    tests/test_diffusion_sampler.py tests/test_diffusion_scope.py \
    tests/test_diffusion_validate.py
  ```

- [x] Run existing model and chirality regression tests:

  ```bash
  pytest -q tests/test_model_fasta_case.py tests/test_model_renumber.py \
    tests/test_model_strip_heterogens.py tests/test_ffutils_chirality.py
  ```

- [ ] Run homology tests when Phase 5 changes begin:

  ```bash
  pytest -q tests/test_homology_plan.py tests/test_homology.py \
    tests/test_homology_fasta.py
  ```

- [x] Run repository documentation and whitespace checks:

  ```bash
  python scripts/check_agent_docs.py
  git diff --check
  ```

- [ ] Run generated-reference updates and checks after argparse changes:

  ```bash
  python scripts/gen_cli_reference.py
  python scripts/gen_gui_spec.py
  python scripts/gen_cli_reference.py --check
  python scripts/gen_gui_spec.py --check
  ```

- [x] Run a pinned GPU smoke test on the reviewed 8CZ8 five-residue internal gap.
- [x] Run a same-seed deterministic repeat on the pinned GPU environment; the raw PDB was byte-identical and coordinate RMSD was `0.0 Å`.
- [x] Run the initial representative RFdiffusion subset and MODELLER comparator: repeated post-refinement seed-7 workspaces in each reviewed 5-/10-residue case are byte-identical, pass every hard gate, and pass the measured MODELLER comparison. Broader strata remain pending.

## First Implementation Iteration Exit Criteria

- [x] Complete all locally implementable Phase 0 and Phase 1 work; externally pinned GPU adapters remain later phases.
- [x] Keep the public CLI unchanged.
- [x] Keep the MODELLER default path unchanged.
- [x] Add no ML dependencies to the core environment.
- [x] Provide a versioned backend-neutral contract.
- [x] Provide deterministic mask construction.
- [x] Provide weighted Kabsch and reinjection primitives.
- [x] Provide a deterministic fake external runner.
- [x] Provide independent validation and atomic no-replace publication.
- [x] Provide a separate diffusion provenance manifest.
- [x] Provide a benchmark manifest with predeclared thresholds.
- [x] Provide CPU-testable adapter preflight and sampler-conformance decisions without claiming a real engine adapter.
- [x] Pass the CPU/Linux macOS-compatible core suite before beginning the RFdiffusion GPU adapter; an actual Apple Silicon run remains a separate platform check.

## Remaining External Gates

- [ ] Phase 2 remains blocked on immutable service-image identity, final checkpoint redistribution review, measured interface evidence, and insertion-code benchmark evidence. Repeatable passing refined candidates and MODELLER comparisons are complete for 8CZ8, independent 8B01 regular loops, and target-chain-only 7X35 glycine/proline-rich difficult loops; the retained-partner 7X35 interface case is declared and locally validated but not yet run on the pinned A100 environment.
- [ ] Phase 3 remains blocked until a pinned all-atom sampler demonstrates externally controllable mutable state, stable atom identity, and exact fixed-coordinate overwrite at every denoising step.
- [ ] Phases 4-6 remain blocked on the Phase 2/3 evidence and intentionally make no public CLI, homology, or production-support claim.
- [ ] Keep status `proposed` and MODELLER/PDBFixer as supported production baselines until those gates pass.
